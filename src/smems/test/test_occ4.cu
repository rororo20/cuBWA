#include "FMI_search.h"
#include "gtest/gtest.h"
#include <map>

#define THREADS_PER_BLOCK 32
void generate_one_hot_mask(unsigned short bwt_mask[][4])
{
    uint64_t *one_hot_mask_array = (uint64_t *)_mm_malloc(64 * sizeof(uint64_t), 64);
    one_hot_mask_array[0] = 0;
    uint64_t base = 0x8000000000000000L;
    one_hot_mask_array[1] = base;
    int64_t i = 0;
    for (i = 2; i < 64; i++) {
        one_hot_mask_array[i] = (one_hot_mask_array[i - 1] >> 1) | base;
    }
    for (i = 0; i < 64; i++) {
        uint64_t offset = one_hot_mask_array[i];
        bwt_mask[i][0] = bwt_mask[i][1] = bwt_mask[i][2] = bwt_mask[i][3] = 0;
        for (int j = 0; j < 16; j++) {
            bwt_mask[i][0] = (bwt_mask[i][0] << 1) | (offset >> 63 & 0X1L);
            bwt_mask[i][1] = (bwt_mask[i][1] << 1) | (offset >> 62 & 0x1L);
            bwt_mask[i][2] = (bwt_mask[i][2] << 1) | (offset >> 61 & 0x1L);
            bwt_mask[i][3] = (bwt_mask[i][3] << 1) | (offset >> 60 & 0x1L);
            offset = offset << 4;
        }
    }
}

void generate_occ_cpp(CP_OCC &cpo, char base[65], int64_t cp_count[4])
{
    uint8_t enc_bases[64];
    cpo.cp_count[0] = cp_count[0];
    cpo.cp_count[1] = cp_count[1];
    cpo.cp_count[2] = cp_count[2];
    cpo.cp_count[3] = cp_count[3];

    for (int i = 0; i < 64; i++) {
        switch (base[i]) {
            case 'A':
                /* code */
                enc_bases[i] = 0;
                break;
            case 'C':
                enc_bases[i] = 1;
                break;
            case 'G':
                enc_bases[i] = 2;
                break;
            case 'T':
                enc_bases[i] = 3;
                break;
            default:
                enc_bases[i] = 4;
                break;
        }
    }
    cpo.one_hot_bwt_str[0] = 0;
    cpo.one_hot_bwt_str[1] = 0;
    cpo.one_hot_bwt_str[2] = 0;
    cpo.one_hot_bwt_str[3] = 0;
    for (int i = 0; i < 64; i++) {
        cpo.one_hot_bwt_str[0] = cpo.one_hot_bwt_str[0] << 1;
        cpo.one_hot_bwt_str[1] = cpo.one_hot_bwt_str[1] << 1;
        cpo.one_hot_bwt_str[2] = cpo.one_hot_bwt_str[2] << 1;
        cpo.one_hot_bwt_str[3] = cpo.one_hot_bwt_str[3] << 1;
        uint8_t c = enc_bases[i];
        if (c < 4) {
            cpo.one_hot_bwt_str[c] += 1;
        }
    }
}
class BackwardTest : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        generate_one_hot_mask(bwt_mask);
        size_t array_size = sizeof(bwt_mask);
        cudaMalloc(&bwt_mask_device, array_size);
        cudaMemcpy(bwt_mask_device, bwt_mask, array_size, cudaMemcpyHostToDevice);
    }
    static void TearDownTestSuite()
    {
        cudaFree(bwt_mask_device);
    }
    void SetUp() {}
    void TearDown() {}
    static unsigned short bwt_mask[64][4];
    static unsigned short *bwt_mask_device;
    // static CP_OCC *cpos;
    // static CP_OCC *cpos_device;
    // static int cpo_size;
};

unsigned short BackwardTest::bwt_mask[64][4];
unsigned short *BackwardTest::bwt_mask_device = NULL;

TEST_F(BackwardTest, testcase1)
{
    EXPECT_EQ(bwt_mask[0][0], 0X0000);
    EXPECT_EQ(bwt_mask[0][1], 0X0000);
    EXPECT_EQ(bwt_mask[0][2], 0X0000);
    EXPECT_EQ(bwt_mask[0][3], 0X0000);

    EXPECT_EQ(bwt_mask[1][0], 0X8000);
    EXPECT_EQ(bwt_mask[1][1], 0X0000);
    EXPECT_EQ(bwt_mask[1][2], 0X0000);
    EXPECT_EQ(bwt_mask[1][3], 0X0000);

    EXPECT_EQ(bwt_mask[2][0], 0X8000);
    EXPECT_EQ(bwt_mask[2][1], 0X8000);
    EXPECT_EQ(bwt_mask[2][2], 0X0000);
    EXPECT_EQ(bwt_mask[2][3], 0X0000);

    EXPECT_EQ(bwt_mask[3][0], 0X8000);
    EXPECT_EQ(bwt_mask[3][1], 0X8000);
    EXPECT_EQ(bwt_mask[3][2], 0X8000);
    EXPECT_EQ(bwt_mask[3][3], 0X0000);

    EXPECT_EQ(bwt_mask[4][0], 0X8000);
    EXPECT_EQ(bwt_mask[4][1], 0X8000);
    EXPECT_EQ(bwt_mask[4][2], 0X8000);
    EXPECT_EQ(bwt_mask[4][3], 0X8000);

    EXPECT_EQ(bwt_mask[5][0], 0XC000);
    EXPECT_EQ(bwt_mask[5][1], 0X8000);
    EXPECT_EQ(bwt_mask[5][2], 0X8000);
    EXPECT_EQ(bwt_mask[5][3], 0X8000);

    EXPECT_EQ(bwt_mask[6][0], 0XC000);
    EXPECT_EQ(bwt_mask[6][1], 0XC000);
    EXPECT_EQ(bwt_mask[6][2], 0X8000);
    EXPECT_EQ(bwt_mask[6][3], 0X8000);

    EXPECT_EQ(bwt_mask[7][0], 0XC000);
    EXPECT_EQ(bwt_mask[7][1], 0XC000);
    EXPECT_EQ(bwt_mask[7][2], 0XC000);
    EXPECT_EQ(bwt_mask[7][3], 0X8000);

    EXPECT_EQ(bwt_mask[8][0], 0XC000);
    EXPECT_EQ(bwt_mask[8][1], 0XC000);
    EXPECT_EQ(bwt_mask[8][2], 0XC000);
    EXPECT_EQ(bwt_mask[8][3], 0XC000);

    EXPECT_EQ(bwt_mask[9][0], 0XE000);
    EXPECT_EQ(bwt_mask[9][1], 0XC000);
    EXPECT_EQ(bwt_mask[9][2], 0XC000);
    EXPECT_EQ(bwt_mask[9][3], 0XC000);

    EXPECT_EQ(bwt_mask[10][0], 0XE000);
    EXPECT_EQ(bwt_mask[10][1], 0XE000);
    EXPECT_EQ(bwt_mask[10][2], 0XC000);
    EXPECT_EQ(bwt_mask[10][3], 0XC000);

    EXPECT_EQ(bwt_mask[11][0], 0XE000);
    EXPECT_EQ(bwt_mask[11][1], 0XE000);
    EXPECT_EQ(bwt_mask[11][2], 0XE000);
    EXPECT_EQ(bwt_mask[11][3], 0XC000);

    EXPECT_EQ(bwt_mask[12][0], 0XE000);
    EXPECT_EQ(bwt_mask[12][1], 0XE000);
    EXPECT_EQ(bwt_mask[12][2], 0XE000);
    EXPECT_EQ(bwt_mask[12][3], 0XE000);

    EXPECT_EQ(bwt_mask[13][0], 0XF000);
    EXPECT_EQ(bwt_mask[13][1], 0XE000);
    EXPECT_EQ(bwt_mask[13][2], 0XE000);
    EXPECT_EQ(bwt_mask[13][3], 0XE000);

    EXPECT_EQ(bwt_mask[14][0], 0XF000);
    EXPECT_EQ(bwt_mask[14][1], 0XF000);
    EXPECT_EQ(bwt_mask[14][2], 0XE000);
    EXPECT_EQ(bwt_mask[14][3], 0XE000);

    EXPECT_EQ(bwt_mask[15][0], 0XF000);
    EXPECT_EQ(bwt_mask[15][1], 0XF000);
    EXPECT_EQ(bwt_mask[15][2], 0XF000);
    EXPECT_EQ(bwt_mask[15][3], 0XE000);

    EXPECT_EQ(bwt_mask[16][0], 0XF000);
    EXPECT_EQ(bwt_mask[16][1], 0XF000);
    EXPECT_EQ(bwt_mask[16][2], 0XF000);
    EXPECT_EQ(bwt_mask[16][3], 0XF000);

    EXPECT_EQ(bwt_mask[60][0], 0XFFFE);
    EXPECT_EQ(bwt_mask[60][1], 0XFFFE);
    EXPECT_EQ(bwt_mask[60][2], 0XFFFE);
    EXPECT_EQ(bwt_mask[60][3], 0XFFFE);

    EXPECT_EQ(bwt_mask[61][0], 0XFFFF);
    EXPECT_EQ(bwt_mask[61][1], 0XFFFE);
    EXPECT_EQ(bwt_mask[61][2], 0XFFFE);
    EXPECT_EQ(bwt_mask[61][3], 0XFFFE);

    EXPECT_EQ(bwt_mask[62][0], 0XFFFF);
    EXPECT_EQ(bwt_mask[62][1], 0XFFFF);
    EXPECT_EQ(bwt_mask[62][2], 0XFFFE);
    EXPECT_EQ(bwt_mask[62][3], 0XFFFE);

    EXPECT_EQ(bwt_mask[63][0], 0XFFFF);
    EXPECT_EQ(bwt_mask[63][1], 0XFFFF);
    EXPECT_EQ(bwt_mask[63][2], 0XFFFF);
    EXPECT_EQ(bwt_mask[63][3], 0XFFFE);
}
TEST_F(BackwardTest, testcase2)
{
    char arr[65] = "ACGTACACGGTTTTTAAAAGGGGCCCCTTTTACGTACGTTGCATGCACGTATTACAACGATTAC";
    /*result*/
    std::map<int, std::vector<int>> result = {
        {0, {1, 0, 0, 0}},      {1, {1, 1, 0, 0}},      {2, {1, 1, 1, 0}},      {3, {1, 1, 1, 1}},      {4, {2, 1, 1, 1}},
        {5, {2, 2, 1, 1}},      {6, {3, 2, 1, 1}},      {7, {3, 3, 1, 1}},      {8, {3, 3, 2, 1}},      {9, {3, 3, 3, 1}},
        {10, {3, 3, 3, 2}},     {11, {3, 3, 3, 3}},     {12, {3, 3, 3, 4}},     {13, {3, 3, 3, 5}},     {14, {3, 3, 3, 6}},
        {15, {4, 3, 3, 6}},     {16, {5, 3, 3, 6}},     {17, {6, 3, 3, 6}},     {18, {7, 3, 3, 6}},     {19, {7, 3, 4, 6}},
        {20, {7, 3, 5, 6}},     {21, {7, 3, 6, 6}},     {22, {7, 3, 7, 6}},     {23, {7, 4, 7, 6}},     {24, {7, 5, 7, 6}},
        {25, {7, 6, 7, 6}},     {26, {7, 7, 7, 6}},     {27, {7, 7, 7, 7}},     {28, {7, 7, 7, 8}},     {29, {7, 7, 7, 9}},
        {30, {7, 7, 7, 10}},    {31, {8, 7, 7, 10}},    {32, {8, 8, 7, 10}},    {33, {8, 8, 8, 10}},    {34, {8, 8, 8, 11}},
        {35, {9, 8, 8, 11}},    {36, {9, 9, 8, 11}},    {37, {9, 9, 9, 11}},    {38, {9, 9, 9, 12}},    {39, {9, 9, 9, 13}},
        {40, {9, 9, 10, 13}},   {41, {9, 10, 10, 13}},  {42, {10, 10, 10, 13}}, {43, {10, 10, 10, 14}}, {44, {10, 10, 11, 14}},
        {45, {10, 11, 11, 14}}, {46, {11, 11, 11, 14}}, {47, {11, 12, 11, 14}}, {48, {11, 12, 12, 14}}, {49, {11, 12, 12, 15}},
        {50, {12, 12, 12, 15}}, {51, {12, 12, 12, 16}}, {52, {12, 12, 12, 17}}, {53, {13, 12, 12, 17}}, {54, {13, 13, 12, 17}},
        {55, {14, 13, 12, 17}}, {56, {15, 13, 12, 17}}, {57, {15, 14, 12, 17}}, {58, {15, 14, 13, 17}}, {59, {16, 14, 13, 17}},
        {60, {16, 14, 13, 18}}, {61, {16, 14, 13, 19}}, {62, {17, 14, 13, 19}}, {63, {17, 15, 13, 19}}};

    CP_OCC *cp = (CP_OCC *)malloc(1024 * sizeof(CP_OCC));
    CP_OCC *cp_device = NULL;
    int64_t count1[] = {0, 0, 0, 0};
    int base_count = 4;
    int N = CP_BLOCK_SIZE * 4;
    generate_occ_cpp(cp[0], arr, count1);
    generate_occ_cpp(cp[1], arr, count1);

    cudaMalloc(&cp_device, 1024 * sizeof(CP_OCC));
    cudaMemcpy(cp_device, cp, 1024 * sizeof(CP_OCC), cudaMemcpyHostToDevice);

    SMEM_CUDA *smems = (SMEM_CUDA *)malloc(sizeof(SMEM_CUDA) * N);
    SMEM_CUDA *smems_device = NULL;

    for (int i = 0; i < N; i++) {
        smems[i].k = i;
        smems[i].l = 0;
        smems[i].s = 0;
    }
    cudaMalloc(&smems_device, N * sizeof(SMEM_CUDA));
    cudaMemcpy(smems_device, smems, N * sizeof(SMEM_CUDA), cudaMemcpyHostToDevice);

    uint8_t *bases = (uint8_t *)malloc(sizeof(uint8_t) * N);
    uint8_t *bases_deviece = NULL;

    for (int i = 0; i < CP_BLOCK_SIZE; i++) {
        for (int j = 0; j < base_count; j++) {
            bases[i * base_count + j] = j;
        }
    }
    cudaMalloc(&bases_deviece, N * sizeof(uint8_t));
    cudaMemcpy(bases_deviece, bases, N * sizeof(uint8_t), cudaMemcpyHostToDevice);

    dim3 block(THREADS_PER_BLOCK);
    dim3 grid(1024 * 1024);
    getOCC4Back<<<grid, block>>>(cp_device, smems_device, bwt_mask_device, bases_deviece, N, 1);

    cudaMemcpy(smems, smems_device, N * sizeof(SMEM_CUDA), cudaMemcpyDeviceToHost);

    for (int i = 1; CP_BLOCK_SIZE; i++) {
        for (int j = 0; j < base_count; j++) {
            EXPECT_EQ(smems[i * base_count + j].k, result[i][j]);
            EXPECT_EQ(smems[i * base_count + j].s, 0);
        }
    }
    cudaFree(smems_device);
    cudaFree(cp_device);
    cudaFree(bases_deviece);
    free(smems);
    free(cp);
    free(bases);
}
