
#ifndef __SMEM_H
#define __SMEM_H
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include "FMI_search.h"
#include <curand_kernel.h>
#include <memory_resource>
typedef struct smem_struct_cuda
{
#ifdef DEBUG
    uint64_t info;  // for debug
#endif
    int64_t k, l, s;
} SMEM_CUDA;

enum class StepStatus { first_pass, second_pass, bwt_seed_strategy };

enum class SearchDirectionStatus { backward, foward, only_foward };

struct SMEMS_STATUS
{
    StepStatus step_status;
    SearchDirectionStatus direct_status;
    int rid;
    int anchor;
    int seq_offset;
    int rightmost;
    // for
    int prev_offset;
    int currr_offset;
    bool has_optimal_smems_occurred;
};

struct FirstPassSmems
{
    int rid;
    int rightmost;
};
class SMEMSerach
{
public:
    SMEMSerach(FMI_search *fm, int smems_size, int reads_size)
        : fmi(fm)
        , batch_reads_size(reads_size)
        , batch_smems_size(smems_size)
        , result_num(0)
        , max_result_num(coefficient_result * batch_reads_size)
        , prev_num(0)
        , max_pre_num(coefficient_prev * batch_reads_size)
    {
        // fmi interface get cp_occ length
        int cp_occ_size = fmi->get_cp_occ_size();
        // copy cp_occ to GPU
        cudaMalloc(&cp_occ, cp_occ_size * sizeof(CP_OCC));
        cudaMemcpy(cp_occ, fmi->get_cp_occ(), cp_occ_size * sizeof(CP_OCC), cudaMemcpyHostToDevice);
        // fmi interface get bwt mask array
        unsigned short *bwt_mask = fmi->get_bwt_mask();
        cudaMalloc(&bwt_mask_device, 64 * 4 * sizeof(unsigned short));
        cudaMemcpy(bwt_mask_device, bwt_mask, 64 * 4 * sizeof(unsigned short), cudaMemcpyHostToDevice);
        // copy cp_occ to GPU
        host_smems = (SMEM_CUDA *)malloc(sizeof(SMEM_CUDA) * batch_smems_size);
        cudaMalloc(&device_smems, sizeof(device_smems) * batch_smems_size);
        // Bases
        host_bases = (uint8_t *)malloc(sizeof(uint8_t) * batch_smems_size);
        cudaMalloc(&device_bases, sizeof(uint8_t) * batch_smems_size);

        result = (SMEM *)malloc(sizeof(SMEM) * max_result_num);
        prev = (SMEM *)malloc(sizeof(SMEM) * max_pre_num);
        status = (SMEMS_STATUS *)malloc(sizeof(SMEMS_STATUS) * batch_smems_size);

        int *running_idx = (int *)malloc(sizeof(int) * batch_smems_size);
        int *idle_idx = (int *)malloc(sizeof(int) * batch_smems_size);
        int *prev_num_ = (int *)malloc(sizeof(int) * batch_smems_size);
    }

    SMEM *collect_smem(const bseq1_t *seq, int nseq, int32_t min_interval);

    ~SMEMSerach()
    {
        cudaFree(cp_occ);
        cudaFree(bwt_mask_device);
        free(host_smems);
        cudaFree(device_smems);
        free(host_bases);
        cudaFree(device_bases);
        free(result);
        free(prev);
        free(status);

        free(running_idx);
        free(idle_idx);
        free(prev_num_);
    }

private:
    static constexpr int thread_per_block = 32;
    static constexpr int block_number = 1024 * 1024;
    static constexpr int coefficient_result = 64;
    static constexpr int coefficient_prev = 150;

    void backward(int process_number);

    int get_maxlength(const bseq1_t *seq, int nseq);

    void initialze_batch_smems_stauts(int &current_seq_id, int &running_reads, const bseq1_t *seq, const int n_seq, const int max_length);

    void non_bwt_seed_foward(const bseq1_t *seq, int &new_running_idx, SMEMS_STATUS *curr, const int smems_idx, const int max_length,
                             const int min_interval);

    void non_bwt_seed_backward(const bseq1_t *seq, int &new_running_idx, int &new_idle_idx, SMEMS_STATUS *curr, const int smems_idx,
                               const int max_length, const int min_interval);

    void bwt_seed_foward(const bseq1_t *seq, int &new_running_idx, int &new_idle_idx, SMEMS_STATUS *curr, const int smems_idx,
                         const int max_length, const int min_interval);

    void prepre_batch_status(const bseq1_t *seq, int &current_seq_id, int &running_reads, int &new_running_idx, int &new_idle_idx,
                             int nseq);

    FMI_search *fmi;
    //  Device Index
    CP_OCC *cp_occ;
    unsigned short *bwt_mask_device;
    // HOST count register for FM-Index
    int64_t count[5];  // TODO: init

    SMEM_CUDA *host_smems;
    SMEM_CUDA *device_smems;

    int64_t sentinel_index;  // TODO: init
    uint8_t *host_bases;
    uint8_t *device_bases;
    int batch_smems_size;
    int batch_reads_size;
    FirstPassSmems *first_result;
    int first_result_num;
    SMEM *result;
    int result_num;
    int max_result_num;
    SMEM *prev;
    int prev_num;
    int max_pre_num;
    SMEMS_STATUS *status;
    int *running_idx;
    int *idle_idx;
    int *prev_num_;
};

__global__ void getOCC4Back(CP_OCC *cp_occ, SMEM_CUDA *smems, unsigned short *bwt_mask, uint8_t *bases, int size, int64_t sentinel_index);

__device__ __host__ uint8_t countSetBits_loop(unsigned short n);

__device__ __host__ uint8_t countSetBits_v1(unsigned short n);

__device__ __host__ uint8_t countSetBits_v2(unsigned short n);
#endif