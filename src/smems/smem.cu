#include "smem.h"
#include "FMI_search.h"

void SMEM_SEARCH::backward(int process_number)
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
