cmake_minimum_required(VERSION 3.10)
project(cuBWA LANGUAGES C CXX CUDA)
set(CMAKE_C_FLAGS "-g -O2 -Wall")
set(CMAKE_CXX_FLAGS "-g -O2 -Wall")
find_package(CUDAToolkit REQUIRED)


include_directories(
    ext/safestringlib/include
)

set(MEM2_FM_INDEX_SOURCE
    memcpy_bwamem.cpp
    FMI_search.cu
    read_index_ele.cpp
    utils.cpp
    bntseq.cpp
)

set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

add_subdirectory(ext/safestringlib)

add_executable(
    cuBWA
    main.cc
    ${MEM2_FM_INDEX_SOURCE}
)

target_include_directories(cuBWA PUBLIC ext/safestringlib/include)


link_libraries(cuBWA
    ext/safestringlib
)

set_target_properties(cuBWA PROPERTIES CUDA_ARCHITECTURES 70)

target_link_libraries(cuBWA 
    PRIVATE 
    safestring_shared
    m 
    z
    CUDA::cudart
)