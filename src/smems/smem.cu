#include "smem.h"
#include "FMI_search.h"

SMEM *SMEMSerach::collect_smem(const bseq1_t *seq, int nseq, int32_t min_interval)
{
    /* init  SMEMS status*/
    int current_seq_id = 0;
    int min_seq_id = std::min(batch_smems_size, nseq);
    int running_reads = 0;
    int max_length = -1;

    for (int i = 0; i < nseq; i++) {
        if (seq[i].l_seq > max_length) {
            max_length = seq[i].l_seq;
        }
    }
    if (max_length * nseq > max_pre_num) {
        prev = (SMEM *)realloc(prev, max_length * sizeof(SMEM));
    }
    // TODO: resource reside in object
    int *running_idx = (int *)malloc(sizeof(int) * batch_smems_size);
    int *idle_idx = (int *)malloc(sizeof(int) * batch_smems_size);
    int *start = (int *)malloc(sizeof(int) * batch_smems_size);
    int *end = (int *)malloc(sizeof(int) * batch_smems_size);
    int *s = (int *)malloc(sizeof(int) * batch_smems_size);
    int *prev_num_ = (int *)malloc(sizeof(int) * batch_smems_size);

    while (current_seq_id < min_seq_id) {
        status[current_seq_id].rid = current_seq_id;
        status[current_seq_id].step_status = StepStatus::first_pass;
        status[current_seq_id].direct_status = SearchDirectionStatus::foward;
        status[current_seq_id].seq_offset = 1;
        current_seq_id++;
        running_reads++;
    }
    /* prepare smems array */
    for (int i = 0; i < running_reads; i++) {
        host_bases[i] = 3 - seq[i].seq[1];
        start[i] = 0;
        end[i] = 0;
        host_smems[i].k = count[3 - seq[i].seq[0]];
        host_smems[i].l = count[seq[i].seq[0]];
        host_smems[i].s = count[seq[i].seq[0] + 1] - count[seq[i].seq[0]];
        s[i] = host_smems[i].l - host_smems[i].k;
        running_idx[i] = i;
        SMEM *tmp = prev + i * max_length;
        tmp->m = tmp->n = 0;
        tmp->rid = i;
        tmp->s = s[i];
        prev_num_[i] = 0;
    }
    int new_idle_idx = 0;
    do {
        // TODO: less  use cpu method
        backward(running_reads);
        // update status
        int new_running_idx = 0;
        for (int i = 0; i < running_reads; i++) {
            SMEMS_STATUS *curr = status + running_idx[i];
            // update host MEMS  and Modify idle_idx
            if (curr->step_status != StepStatus::bwt_seed_strategy) {
                if (curr->direct_status == SearchDirectionStatus::foward) {
                    // todo: l_seq
                    int mapping_prev_offset = i * max_length;
                    //  before backward may save
                    prev_num_[i] += host_smems[i].s != prev[mapping_prev_offset + prev_num_[i]].s;

                    prev[mapping_prev_offset + prev_num_[i]].s = host_smems[i].s;
                    // SWAP
                    prev[mapping_prev_offset + prev_num_[i]].k = host_smems[i].l;
                    prev[mapping_prev_offset + prev_num_[i]].l = host_smems[i].k;
                    prev[mapping_prev_offset + prev_num_[i]].n = curr->seq_offset;
                    if (host_smems[i].s < min_interval || curr->seq_offset == seq[curr->rid].l_seq - 1) {  // switch to BackWard
                        // Last interval
                        if (curr->seq_offset != seq[curr->rid].l_seq - 1)
                            prev_num_[i] += (prev[mapping_prev_offset + prev_num_[i]].s >= min_interval);
                        // SORT ...
                        for (int p = 0; p < (prev_num_[i] / 2); p++) {
                            SMEM temp = prev[mapping_prev_offset + p];
                            prev[mapping_prev_offset + p] = prev[mapping_prev_offset + prev_num_[i] - p - 1];
                            prev[mapping_prev_offset + prev_num_[i] - p - 1] = temp;
                        }
                        curr->rightmost = curr->seq_offset;
                        if (curr->anchor == 0 || prev_num_[i] == 0) {  // continue  foward search from next anchor's

                            curr->anchor = curr->seq_offset + 1;
                            host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->anchor + 1];
                            host_smems[new_running_idx].k = count[3 - seq[curr->rid].seq[curr->anchor]];
                            host_smems[new_running_idx].l = count[seq[curr->rid].seq[curr->anchor]];
                            host_smems[new_running_idx].s =
                                count[seq[curr->rid].seq[curr->anchor] + 1] - count[seq[curr->rid].seq[curr->anchor]];
                            curr->seq_offset = curr->anchor + 2;
                            running_idx[new_running_idx] = running_idx[i];
                            new_running_idx++;
                            curr->direct_status = SearchDirectionStatus::foward;
                        }
                        else {
                            // starting backward search
                            curr->seq_offset = curr->anchor - 1;
                            host_bases[new_running_idx] = seq[curr->rid].seq[curr->seq_offset];
                            running_idx[new_running_idx] = running_idx[i];
                            host_smems[new_running_idx].k = prev[mapping_prev_offset].k;
                            host_smems[new_running_idx].l = prev[mapping_prev_offset].l;
                            host_smems[new_running_idx].s = prev[mapping_prev_offset].s;
                            new_running_idx++;
                            curr->currr_offset = 0;
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
                        running_idx[new_running_idx] = running_idx[i];
                        new_running_idx++;
                    }
                }
                else {
                    if (curr->prev_offset == prev_num_[i] - 1) {
                        if (curr->currr_offset == 0) {  // backward can't extend , switch new step
                            if (curr->rightmost == seq[curr->rid].l_seq - 1) {
                                if (curr->step_status == StepStatus::first_pass) {
                                    idle_idx[new_idle_idx] = running_idx[i];
                                    new_idle_idx++;
                                }
                                else {
                                    // swith bwt seeds
                                    curr->step_status = StepStatus::bwt_seed_strategy;
                                }
                            }
                            else {  // Switch forward
                                curr->anchor = curr->rightmost;
                                host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->anchor + 1];
                                host_smems[new_running_idx].k = count[3 - seq[curr->rid].seq[curr->anchor]];
                                host_smems[new_running_idx].l = count[seq[curr->rid].seq[curr->anchor]];
                                host_smems[new_running_idx].s =
                                    count[seq[curr->rid].seq[curr->anchor] + 1] - count[seq[curr->rid].seq[curr->anchor]];
                                curr->seq_offset = curr->rightmost + 1;
                                running_idx[new_running_idx] = running_idx[i];
                                new_running_idx++;
                                ;
                            }
                        }
                        else {  // continue backward
                            curr->seq_offset--;
                            host_bases[new_running_idx] = seq[curr->rid].seq[curr->seq_offset];
                            host_smems[new_running_idx].k = prev[i * max_length].k;
                            host_smems[new_running_idx].l = prev[i * max_length].l;
                            host_smems[new_running_idx].s = prev[i * max_length].s;
                            prev_num_[i] = curr->currr_offset;
                            curr->currr_offset = 0;
                            curr->prev_offset = 0;
                            running_idx[new_running_idx] = running_idx[i];
                            new_running_idx++;
                        }
                    }
                    else {
                        if (prev[i * max_length + curr->prev_offset].s != host_smems[i].s) {
                            if (!curr->has_optimal_smems_occurred && host_smems[i].s < min_interval &&

                                curr->rightmost - curr->seq_offset > 0) {
                                // push to results
                                curr->has_optimal_smems_occurred = true;
                                result[result_num].k = host_smems[i].k;
                                result[result_num].l = host_smems[i].l;
                                result[result_num].s = host_smems[i].s;
                                result[result_num].rid = curr->rid;
                                result[result_num].m = curr->seq_offset;
                                result[result_num].n = prev[i * max_length + curr->prev_offset].n;
                                // push rightmost
                                first_result[first_result_num].rid = curr->rid;
                                first_result[first_result_num].rightmost = curr->rightmost;
                                first_result_num++;
                                result_num++;
                            }

                            if (host_smems[i].s >= min_interval) {
                                prev[i * max_length + curr->currr_offset].k = host_smems[i].k;
                                prev[i * max_length + curr->currr_offset].l = host_smems[i].l;
                                prev[i * max_length + curr->currr_offset].s = host_smems[i].s;
                                curr->currr_offset++;
                            }
                        }
                        curr->prev_offset++;
                        host_bases[new_running_idx] = seq[curr->rid].seq[curr->seq_offset];
                        host_smems[new_running_idx].k = prev[i * max_length + curr->prev_offset].k;
                        host_smems[new_running_idx].l = prev[i * max_length + curr->prev_offset].l;
                        host_smems[new_running_idx].s = prev[i * max_length + curr->prev_offset].s;
                        running_idx[new_running_idx] = running_idx[i];
                        new_running_idx++;
                    }
                }
            }
            else {
                curr->seq_offset++;
                if (host_smems[new_running_idx].s < min_interval) {
                    ;
                }
                host_bases[new_running_idx] = 3 - seq[curr->rid].seq[curr->seq_offset];
                host_smems[new_running_idx].k = host_smems[i].l;
                host_smems[new_running_idx].l = host_smems[i].k;
                host_smems[new_running_idx].s = host_smems[i].s;
                running_idx[new_running_idx] = running_idx[i];
            }
        }
        //  prepare new reads's array
        while (new_running_idx < batch_smems_size && current_seq_id < nseq) {
            int new_idx = idle_idx[new_idle_idx];
            if (first_result_num > 0) {
                int rid = first_result[first_result_num].rid;
                int anchor = first_result[first_result_num].rightmost;
                status[new_idx].rid = rid;
                status[new_idx].step_status = StepStatus::second_pass;
                status[new_idx].direct_status = SearchDirectionStatus::foward;
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
                status[new_idx].direct_status = SearchDirectionStatus::foward;
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
    } while (running_reads > 0);
    free(running_idx);
    free(start);
    free(end);
    free(s);
    return NULL;
}
/* TODO: Stream*/
void SMEMSerach::backward(int process_number)
{
    if (process_number == 0) return;
    cudaMemcpy(device_bases, host_bases, process_number * sizeof(uint8_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_smems, host_smems, process_number * sizeof(SMEM_CUDA), cudaMemcpyHostToDevice);
    dim3 block(thread_per_block);
    dim3 grid(block_number);
    getOCC4Back<<<grid, block>>>(cp_occ, device_smems, bwt_mask_device, device_bases, process_number, sentinel_index);
    cudaMemcpy(host_smems, device_smems, process_number * sizeof(SMEM_CUDA), cudaMemcpyDeviceToHost);

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
