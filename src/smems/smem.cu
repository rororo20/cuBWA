#include "smem.h"
#include "FMI_search.h"

#define CUDA_CHECK(call)                                                                                 \
    do {                                                                                                 \
        cudaError_t error = call;                                                                        \
        if (error != cudaSuccess) {                                                                      \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
            exit(EXIT_FAILURE);                                                                          \
        }                                                                                                \
    } while (0)

// Get the maximum length among all sequences
int SMEMSearch::get_maxlength(const bseq1_t *seq, int nseq)
{
    int max_length = -1;

    for (int i = 0; i < nseq; i++) {
        if (seq[i].l_seq > max_length) {
            max_length = seq[i].l_seq;
        }
    }
    return max_length;
}

// Initialize the first batch of reads status
// Key steps:
// 1. Set initial search direction as forward
// 2. Initialize the first SMEM interval information
// 3. Set initial anchor position
// 4. Handle N bases in sequence
void SMEMSearch::initialize_batch_smems_status(int &current_seq_id, int &running_reads, const bseq1_t *seq, const int n_seq,
                                               const int max_length)
{
    int batch_max_number = std::min(batch_smems_size, n_seq);
    while (current_seq_id < batch_max_number) {
        int next_i = 0;
        while (seq[current_seq_id].seq[next_i] >= 4 || seq[current_seq_id].seq[next_i + 1] >= 4) next_i++;
        uint8_t c = seq[current_seq_id].seq[next_i];
        SMEM *tmp = prev + (current_seq_id * max_length);
        status[current_seq_id].rid = current_seq_id;
        status[current_seq_id].step_status = StepStatus::first_pass;
        status[current_seq_id].direct_status = SearchDirectionStatus::forward;
        status[current_seq_id].seq_offset = 1;
        host_bases[current_seq_id] = 3 - seq[current_seq_id].seq[1];
        host_smems[current_seq_id].k = count[3 - c];
        host_smems[current_seq_id].l = count[c];
        host_smems[current_seq_id].s = count[c + 1] - count[c];
        tmp->m = 0;
        tmp->n = 1;
        tmp->rid = current_seq_id;
        tmp->s = host_smems[current_seq_id].s;
        prev_num_[current_seq_id] = 0;
        current_seq_id++;
        running_reads++;
    }
}

// Forward search without BWT seed strategy
// Main logic:
// 1. Extend current SMEM
// 2. Switch to backward search if interval is too small or reach sequence end
// 3. Find new anchor if current anchor is 0 or no valid interval
// 4. Handle N bases and special cases
void SMEMSearch::non_bwt_seed_forward(const bseq1_t *seq, int &new_running_idx, SMEMS_STATUS *curr, const int smems_idx,
                                      const int max_length, const int min_interval)
{
    int mapping_prev_offset = smems_idx * max_length;
    prev_num_[smems_idx] += host_smems[smems_idx].s != prev[mapping_prev_offset + prev_num_[smems_idx]].s;
    SMEM *current_pre = prev + mapping_prev_offset + prev_num_[smems_idx];

    current_pre->s = host_smems[smems_idx].s;
    current_pre->k = host_smems[smems_idx].l;
    current_pre->l = host_smems[smems_idx].k;
    current_pre->n = curr->seq_offset;
    current_pre->m = curr->anchor;
    if (host_smems[smems_idx].s < min_interval || curr->seq_offset == seq[curr->rid].l_seq - 1) {  // switch to BackWard
        // TODO: next offset is N base
        // Last interval
        if (curr->seq_offset != seq[curr->rid].l_seq - 1)
            prev_num_[smems_idx] += (prev[mapping_prev_offset + prev_num_[smems_idx]].s >= min_interval);
        // SORT ...
        for (int p = 0; p < (prev_num_[smems_idx] / 2); p++) {
            SMEM temp = prev[mapping_prev_offset + p];
            prev[mapping_prev_offset + p] = prev[mapping_prev_offset + prev_num_[smems_idx] - p - 1];
            prev[mapping_prev_offset + prev_num_[smems_idx] - p - 1] = temp;
        }
        curr->rightmost = curr->seq_offset;
        if (curr->anchor == 0 || prev_num_[smems_idx] == 0) {  // continue  foward search from next anchor's
            while (seq[curr->rid].seq[curr->seq_offset] >= 4 || seq[curr->rid].seq[curr->seq_offset + 1] >= 4) curr->seq_offset++;
            curr->anchor = curr->seq_offset + 1;
            uint8_t c = seq[curr->rid].seq[curr->anchor];
            host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->anchor + 1];
            host_smems[new_running_idx].k = count[3 - c];
            host_smems[new_running_idx].l = count[c];
            host_smems[new_running_idx].s = count[c + 1] - count[c];
            curr->seq_offset = curr->anchor + 1;
            running_idx[new_running_idx] = running_idx[smems_idx];
            new_running_idx++;
            curr->direct_status = SearchDirectionStatus::forward;
        }
        else {
            // starting backward search
            curr->seq_offset = curr->anchor - 1;
            // TODO: consider N base
            host_bases[new_running_idx] = seq[curr->rid].seq[curr->seq_offset];
            running_idx[new_running_idx] = running_idx[smems_idx];
            host_smems[new_running_idx].k = prev[mapping_prev_offset].k;
            host_smems[new_running_idx].l = prev[mapping_prev_offset].l;
            host_smems[new_running_idx].s = prev[mapping_prev_offset].s;
            new_running_idx++;
            curr->curr_offset = 0;
            curr->prev_offset = 0;
            curr->direct_status = SearchDirectionStatus::backward;
        }
    }
    else {
        //  continue Foward.
        curr->seq_offset++;
        host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->seq_offset];
        host_smems[new_running_idx].k = host_smems[mapping_prev_offset].l;
        host_smems[new_running_idx].l = host_smems[mapping_prev_offset].k;
        host_smems[new_running_idx].s = host_smems[mapping_prev_offset].s;
        running_idx[new_running_idx] = running_idx[smems_idx];
        new_running_idx++;
    }
}

// Backward search without BWT seed strategy
// Key points:
// 1. Check if need to switch search phase
// 2. Switch to new phase if backward extension fails
// 3. Save optimal SMEMs when found
// 4. Handle boundary conditions
void SMEMSearch::non_bwt_seed_backward(const bseq1_t *seq, int &new_running_idx, int &new_idle_idx, SMEMS_STATUS *curr, const int smems_idx,
                                       const int max_length, const int min_interval)
{
    if (curr->prev_offset == prev_num_[smems_idx] - 1) {
        if (curr->curr_offset == 0) {  // backward can't extend , switch new step
            if (curr->rightmost == seq[curr->rid].l_seq - 1) {
                if (curr->step_status == StepStatus::first_pass) {
                    idle_idx[new_idle_idx] = running_idx[smems_idx];
                    new_idle_idx++;
                }
                else {
                    // swith bwt seeds
                    curr->step_status = StepStatus::bwt_seed_strategy;
                    int next_i = 0;
                    while (seq[curr->rid].seq[next_i] >= 4 || seq[curr->rid].seq[next_i + 1] >= 4) next_i++;
                    uint8_t c = seq[curr->rid].seq[next_i + 1];
                    host_bases[new_running_idx] = 3 - seq[curr->rid].seq[next_i];
                    host_smems[new_running_idx].k = count[c];
                    host_smems[new_running_idx].l = count[3 - c];
                    host_smems[new_running_idx].s = count[c + 1] - count[c];
                    curr->seq_offset = next_i + 1;
                    new_running_idx++;
                }
            }
            else {  // Switch forward
                curr->anchor = curr->rightmost;
                host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->anchor + 1];
                host_smems[new_running_idx].k = count[3 - seq[curr->rid].seq[curr->anchor]];
                host_smems[new_running_idx].l = count[seq[curr->rid].seq[curr->anchor]];
                host_smems[new_running_idx].s = count[seq[curr->rid].seq[curr->anchor] + 1] - count[seq[curr->rid].seq[curr->anchor]];
                curr->seq_offset = curr->rightmost + 1;
                running_idx[new_running_idx] = running_idx[smems_idx];
                new_running_idx++;
            }
        }
        else {  // continue backward
            curr->seq_offset--;
            host_bases[new_running_idx] = seq[curr->rid].seq[curr->seq_offset];
            host_smems[new_running_idx].k = prev[smems_idx * max_length].k;
            host_smems[new_running_idx].l = prev[smems_idx * max_length].l;
            host_smems[new_running_idx].s = prev[smems_idx * max_length].s;
            prev_num_[smems_idx] = curr->curr_offset;
            curr->curr_offset = 0;
            curr->prev_offset = 0;
            running_idx[new_running_idx] = running_idx[smems_idx];
            new_running_idx++;
        }
    }
    else {
        if (prev[smems_idx * max_length + curr->prev_offset].s != host_smems[smems_idx].s) {
            if (!curr->has_optimal_smems_occurred && host_smems[smems_idx].s < min_interval && curr->rightmost - curr->seq_offset > 0) {
                // push to results
                curr->has_optimal_smems_occurred = true;
                result[result_num].k = host_smems[smems_idx].k;
                result[result_num].l = host_smems[smems_idx].l;
                result[result_num].s = host_smems[smems_idx].s;
                result[result_num].rid = curr->rid;
                result[result_num].m = curr->seq_offset;
                result[result_num].n = prev[smems_idx * max_length + curr->prev_offset].n;
                // push rightmost
                first_result[first_result_num].rid = curr->rid;
                first_result[first_result_num].rightmost = curr->rightmost;
                first_result_num++;
                result_num++;
            }

            if (host_smems[smems_idx].s >= min_interval) {
                prev[smems_idx * max_length + curr->curr_offset].k = host_smems[smems_idx].k;
                prev[smems_idx * max_length + curr->curr_offset].l = host_smems[smems_idx].l;
                prev[smems_idx * max_length + curr->curr_offset].s = host_smems[smems_idx].s;
                curr->curr_offset++;
            }
        }
        curr->prev_offset++;
        host_bases[new_running_idx] = seq[curr->rid].seq[curr->seq_offset];
        host_smems[new_running_idx].k = prev[smems_idx * max_length + curr->prev_offset].k;
        host_smems[new_running_idx].l = prev[smems_idx * max_length + curr->prev_offset].l;
        host_smems[new_running_idx].s = prev[smems_idx * max_length + curr->prev_offset].s;
        running_idx[new_running_idx] = running_idx[smems_idx];
        new_running_idx++;
    }
}

// Forward search with BWT seed strategy
// Process:
// 1. Extend current SMEM
// 2. Save results and find new anchor if interval is too small
// 3. Continue forward extension otherwise
// 4. Handle special cases and N bases
void SMEMSearch::bwt_seed_forward(const bseq1_t *seq, int &new_running_idx, int &new_idle_idx, SMEMS_STATUS *curr, const int smems_idx,
                                  const int max_length, const int min_interval)
{
    curr->seq_offset++;
    if (host_smems[new_running_idx].s < min_interval) {
        if (prev[smems_idx * max_length].s > 0) {  // store
            result[result_num].k = prev[smems_idx * max_length].k;
            result[result_num].l = prev[smems_idx * max_length].l;
            result[result_num].s = prev[smems_idx * max_length].s;
            result[result_num].rid = curr->rid;
            result[result_num].m = curr->seq_offset - 1;
            result[result_num].n = curr->anchor;
        }
        curr->anchor = curr->seq_offset + 1;
        host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->anchor + 1];
        host_smems[new_running_idx].k = count[3 - seq[curr->rid].seq[curr->anchor]];
        host_smems[new_running_idx].l = count[seq[curr->rid].seq[curr->anchor]];
        host_smems[new_running_idx].s = count[seq[curr->rid].seq[curr->anchor] + 1] - count[seq[curr->rid].seq[curr->anchor]];
        curr->seq_offset = curr->anchor + 1;
        running_idx[new_running_idx] = running_idx[smems_idx];
        new_running_idx++;
    }
    else {
        host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->seq_offset];
        // FIXME: swap  i == new_running_idx
        host_smems[new_running_idx].k = host_smems[smems_idx].l;
        host_smems[new_running_idx].l = host_smems[smems_idx].k;
        host_smems[new_running_idx].s = host_smems[smems_idx].s;
        prev[smems_idx * max_length].k = host_smems[new_running_idx].k;
        prev[smems_idx * max_length].l = host_smems[new_running_idx].l;
        prev[smems_idx * max_length].s = host_smems[new_running_idx].s;
        running_idx[new_running_idx] = running_idx[smems_idx];
    }
    new_running_idx++;
}

// Prepare status for next batch of reads
// Steps:
// 1. Process first pass results
// 2. Initialize new reads
// 3. Update running reads list
// 4. Handle idle indices
void SMEMSearch::prepare_batch_status(const bseq1_t *seq, int &current_seq_id, int &running_reads, int &new_running_idx, int &new_idle_idx,
                                      int nseq)
{
    while (new_running_idx < batch_smems_size && (current_seq_id < nseq || first_result_num > 0)) {
        int new_idx = idle_idx[new_idle_idx];
        if (first_result_num > 0) {
            int rid = first_result[first_result_num].rid;
            int anchor = first_result[first_result_num].rightmost;
            status[new_idx].rid = rid;
            status[new_idx].step_status = StepStatus::second_pass;
            status[new_idx].direct_status = SearchDirectionStatus::forward;
            status[new_idx].seq_offset = anchor + 1;
            host_bases[new_idx] = 3 - seq[rid].seq[anchor + 1];
            host_smems[new_idx].k = count[seq[rid].seq[anchor]];
            host_smems[new_idx].l = count[3 - seq[rid].seq[anchor]];
            host_smems[new_idx].s = count[seq[rid].seq[anchor] + 1] - count[seq[rid].seq[anchor]];
            first_result_num--;
        }
        else {
            status[new_idx].rid = current_seq_id;
            status[new_idx].step_status = StepStatus::first_pass;
            status[new_idx].direct_status = SearchDirectionStatus::forward;
            status[new_idx].seq_offset = 1;
            host_bases[new_idx] = 3 - seq[current_seq_id].seq[1];
            host_smems[new_idx].k = count[seq[current_seq_id].seq[0]];
            host_smems[new_idx].l = count[3 - seq[current_seq_id].seq[0]];
            host_smems[new_idx].s = count[seq[current_seq_id].seq[0] + 1] - count[seq[current_seq_id].seq[0]];
            current_seq_id++;
        }
        running_idx[new_running_idx] = new_idx;
        new_running_idx++;
        new_idle_idx--;
    }
    running_reads = new_running_idx;
}

// Main SMEM collection function
// Workflow:
// 1. Initialize status and memory
// 2. Process all reads until completion
// 3. Choose search strategy based on read status
// 4. Handle results and memory management
SMEM *SMEMSearch::collect_smem(const bseq1_t *seq, int nseq, int32_t min_interval)
{
    /* init  SMEMS status*/
    int current_seq_id = 0;
    int running_reads = 0;
    int max_length = get_maxlength(seq, nseq);

    if (max_length * nseq > max_pre_num) {
        prev = (SMEM *)realloc(prev, max_length * sizeof(SMEM));
        if (prev == nullptr) {
            printf("Error: Failed to reallocate prev array\n");
            exit(1);
        }
        max_pre_num = max_length * nseq;
    }
    initialize_batch_smems_status(current_seq_id, running_reads, seq, nseq, max_length);
    /* prepare smems array */
    int new_idle_idx = 0;
    do {
        // TODO: less  use cpu method
        backward(running_reads);
        int new_running_idx = 0;
        for (int i = 0; i < running_reads; i++) {
            SMEMS_STATUS *curr = status + running_idx[i];
            // update host MEMS  and Modify idle_idx
            if (curr->step_status != StepStatus::bwt_seed_strategy) {
                if (curr->direct_status == SearchDirectionStatus::forward) {
                    non_bwt_seed_forward(seq, new_running_idx, curr, i, max_length, min_interval);
                }
                else {
                    non_bwt_seed_backward(seq, new_running_idx, new_idle_idx, curr, i, max_length, min_interval);
                }
            }
            else {
                bwt_seed_forward(seq, new_running_idx, new_idle_idx, curr, i, max_length, min_interval);
            }
        }
        //  prepare new reads's array
        prepare_batch_status(seq, current_seq_id, running_reads, new_running_idx, new_idle_idx, nseq);
    } while (running_reads > 0);
    return NULL;
}

// GPU implementation of backward search
// Process:
// 1. Copy data to GPU
// 2. Execute GPU kernel computation
// 3. Get results and update count
// 4. Handle CUDA errors
void SMEMSearch::backward(int process_number)
{
    if (process_number == 0) return;
    CUDA_CHECK(cudaMemcpy(device_bases, host_bases, process_number * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(device_smems, host_smems, process_number * sizeof(SMEM_CUDA), cudaMemcpyHostToDevice));
    dim3 block(thread_per_block);
    dim3 grid((process_number + block.x - 1) / block.x);
    getOCC4Back<<<grid, block>>>(cp_occ, device_smems, bwt_mask_device, device_bases, process_number, sentinel_index);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(host_smems, device_smems, process_number * sizeof(SMEM_CUDA), cudaMemcpyDeviceToHost));

    for (int i = 0; i < process_number; i++) {
        host_smems[i].k += count[host_bases[i]];
    }
}
#define IS_K(tid)                (tid < 15)
#define GET_GROUP_THREAD_ID(tid) (tid & 0x3)
#define GET_BASE_PAIR(tid)       ((tid >> 2) & 0x3)

/**
 * @brief backwardExt GPU implement
 * @param cp_occ  checkpoint occ scalar
 * @param bwt_mask  bwt mask array
 * @param bases  multipy base pairs 2bit encode
 * @param size  number of base pairs
 * @param sentinel_index  sentinel index in suffix array
 * @return void
 */
__global__ void getOCC4Back(CP_OCC *cp_occ, SMEM_CUDA *smems, unsigned short *bwt_mask, uint8_t *bases, int size, int64_t sentinel_index)
{
    int base_idx = blockIdx.x;
    int tid = threadIdx.x;
    unsigned short mask = 0;
    unsigned short onehot = 0;
    uint8_t count = 0;
    __shared__ int64_t k[4], l[4], s[4];
    if (base_idx >= size) {
        return;
    }
    uint8_t base = bases[base_idx];
    SMEM_CUDA curr = smems[base_idx];
    if (IS_K(tid)) {
        mask = bwt_mask[((curr.k & CP_MASK) << 2) + GET_GROUP_THREAD_ID(tid)];
        onehot = cp_occ[curr.k >> CP_SHIFT].one_hot_bwt_str[GET_BASE_PAIR(tid)] >> ((3 - GET_GROUP_THREAD_ID(tid)) << 4);
    }

    else {
        mask = bwt_mask[(((curr.k + curr.s) & CP_MASK) << 2) + GET_GROUP_THREAD_ID(tid)];
        onehot = cp_occ[(curr.k + curr.s) >> CP_SHIFT].one_hot_bwt_str[GET_BASE_PAIR(tid)] >> ((3 - GET_GROUP_THREAD_ID(tid)) << 4);
    }

    onehot = onehot & mask;
    /*TODO: bit opt */
    for (int i = 0; i < 16; i++) {
        count += (onehot & 0x1);
        onehot = onehot >> 1;
    }
    unsigned int wrap_mask = (0xFu) << ((tid >> 2) << 2);  //  Generate Wrap Mask for sync add count

    count += __shfl_xor_sync(wrap_mask, count, 1);
    count += __shfl_xor_sync(wrap_mask, count, 2);
    if ((GET_GROUP_THREAD_ID(tid)) == 0) {
        if (IS_K(tid)) {  // update k
            k[GET_BASE_PAIR(tid)] = count + cp_occ[curr.k >> CP_SHIFT].cp_count[GET_BASE_PAIR(tid)];
        }
        else {  // update L
            l[GET_BASE_PAIR(tid)] = count + cp_occ[(curr.k + curr.s) >> CP_SHIFT].cp_count[GET_BASE_PAIR(tid)];
        }
    }
    __syncthreads();
    if (tid < 4) {
        s[tid] = l[tid] - k[tid];
    }
    __syncthreads();
    if (tid == 0) {
        l[3] = curr.l + ((curr.k <= sentinel_index) && ((curr.k + curr.s) > sentinel_index));
        l[2] = l[3] + s[3];
        l[1] = l[2] + s[2];
        l[0] = l[1] + s[1];
        smems[base_idx].k = k[base];  // Need add COUNT
        smems[base_idx].l = l[base];
        smems[base_idx].s = s[base];
    }
    __syncthreads();
}

// Bit counting helper functions
// Three different implementations for performance comparison
__device__ __host__ uint8_t countSetBits_loop(unsigned short n)
{
    uint8_t count = 0;
    for (int i = 0; i < 16; i++) {
        count += (n & 0x1);
        n = n >> 1;
    }
    return count;
}

__device__ __host__ uint8_t countSetBits_v1(unsigned short n)
{
    n = (n & 0x5555) + ((n >> 1) & 0x5555);
    n = (n & 0x3333) + ((n >> 2) & 0x3333);
    n = (n & 0x0f0f) + ((n >> 4) & 0x0f0f);
    n = (n & 0x00ff) + ((n >> 8) & 0x00ff);
    return n & 0x001f;
}

__device__ __host__ uint8_t countSetBits_v2(unsigned short n)
{
    n = (n & 0x5555) + ((n >> 1) & 0x5555);  // 2bit * 8
    n = (n & 0x3333) + ((n >> 2) & 0x3333);  // 4bit * 4
    n = (n & 0x0f0f) + ((n >> 4) & 0x0f0f);  // 8bit * 2
    return (n * 0x0101) >> 8;
}
