cmake_minimum_required(VERSION 3.10)
project(cuBWA)
set(CMAKE_OSX_ARCHITECTURES arm64)

set(CMAKE_C_FLAGS "-g -O2 -Wall")
set(CMAKE_CXX_FLAGS "-g -O2 -Wall")


set()

add_executable(
    cuBWA
    main.cc
)