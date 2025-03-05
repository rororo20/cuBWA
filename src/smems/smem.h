
#ifndef __SMEM_H
#define __SMEM_H
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include "FMI_search.h"
#include <curand_kernel.h>

typedef struct smem_struct_cuda
{
#ifdef DEBUG
    uint64_t info;  // for debug
#endif
    int64_t k, l, s;
} SMEM_CUDA;

class SMEM_SEARCH
{
public:
    SMEM_SEARCH(FMI_search *fm, int smems_size)
        : fmi(fm)
        , batch_smems_size(smems_size)
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
    }

    ~SMEM_SEARCH()
    {
        cudaFree(cp_occ);
        cudaFree(bwt_mask_device);
        free(host_smems);
        cudaFree(device_smems);
        free(host_bases);
        cudaFree(device_bases);
    }

private:
    static constexpr int thread_per_block = 32;
    static constexpr int block_number = 1024 * 1024;

    void backward(int process_number);

    FMI_search *fmi;
    //  Device Index
    CP_OCC *cp_occ;
    unsigned short *bwt_mask_device;
    // HOST count register for FM-Index
    int64_t count[5];

    SMEM_CUDA *host_smems;
    SMEM_CUDA *device_smems;

    int64_t sentinel_index;
    uint8_t *host_bases;
    uint8_t *device_bases;
    int batch_smems_size;
};

__global__ void getOCC4Back(CP_OCC *cp_occ, SMEM_CUDA *smems, unsigned short *bwt_mask, uint8_t *bases, int size, int64_t sentinel_index);

__device__ __host__ uint8_t countSetBits_loop(unsigned short n);

__device__ __host__ uint8_t countSetBits_v1(unsigned short n);

__device__ __host__ uint8_t countSetBits_v2(unsigned short n);
#endif