#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#define DATA_SIZE         (1024 * 1024 * 512)  // 测试512MB数据
#define THREADS_PER_BLOCK 256
#define ITERATIONS        100  // 每个线程执行操作次数

// 生成随机访问模式的核函数
__global__ void generate_random_indices(unsigned int* indices, int size, unsigned long long seed)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= size) return;

    curandState_t state;
    curand_init(seed, tid, 0, &state);
    indices[tid] = curand(&state) % size;
}

// 随机读写测试核函数
__global__ void random_access_kernel(float* data, unsigned int* indices, int size)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= size) return;

#pragma unroll
    for (int i = 0; i < ITERATIONS; ++i) {
        // 随机读取 -> 计算 -> 随机写入
        float value = data[indices[tid]];
        value = value * 2.0f + 1.0f;  // 模拟简单计算
        data[indices[tid]] = value;
    }
}

int main()
{
    float* d_data;
    unsigned int* d_indices;
    cudaEvent_t start, stop;
    float elapsed_time;

    // 分配设备内存
    cudaMalloc(&d_data, DATA_SIZE * sizeof(float));
    cudaMalloc(&d_indices, DATA_SIZE * sizeof(unsigned int));

    // 初始化数据
    cudaMemset(d_data, 0, DATA_SIZE * sizeof(float));

    // 创建CUDA事件
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // 生成随机索引
    dim3 block(THREADS_PER_BLOCK);
    dim3 grid((DATA_SIZE + block.x - 1) / block.x);
    generate_random_indices<<<grid, block>>>(d_indices, DATA_SIZE, time(NULL));
    cudaDeviceSynchronize();

    // 执行测试
    cudaEventRecord(start);
    random_access_kernel<<<grid, block>>>(d_data, d_indices, DATA_SIZE);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time, start, stop);

    // 计算带宽
    size_t bytes_accessed = DATA_SIZE * sizeof(float) * 2 * ITERATIONS;  // 读+写
    double bandwidth = (bytes_accessed / (elapsed_time / 1000.0)) / (1024.0 * 1024.0 * 1024.0);

    printf("Elapsed Time: %.3f ms\n", elapsed_time);
    printf("Effective Bandwidth: %.2f GB/s\n", bandwidth);

    // 清理资源
    cudaFree(d_data);
    cudaFree(d_indices);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return 0;
}
