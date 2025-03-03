set(PACK_NAME smems_test)
# find_package(CUDAToolkit REQUIRED)

set(PACK_SRC test_occ4.cu)
set(PACK_INCLUDE_DIR ..)
# set(PACK_LINK_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PACK_LINK_LIB gtest_main gtest smems z CUDA::cudart)


add_executable(${PACK_NAME} ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
find_package(CUDA  REQUIRED)
include_directories("${CUDA_INCLUDE_DIRS}")
# target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

add_test(NAME ${PACK_NAME} COMMAND ${PACK_NAME})
set_target_properties(${PACK_NAME} PROPERTIES CUDA_ARCHITECTURES 70)