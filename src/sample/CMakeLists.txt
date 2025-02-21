cmake_minimum_required(VERSION 3.4)

project(test_band_test LANGUAGES CXX CUDA)



set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)


add_executable(rand_bandwith_read test_random_read.cu)

set_target_properties(rand_bandwith_read PROPERTIES CUDA_ARCHITECTURES 70)