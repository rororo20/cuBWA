set(PACK_NAME smems)
# find_package(CUDAToolkit REQUIRED)
aux_source_directory(. PACK_SRC)
add_library(${PACK_NAME} SHARED ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ../ext/safestringlib/include)
target_link_libraries(${PACK_NAME}  safestring_shared  z CUDA::cudart)
set_target_properties(${PACK_NAME} PROPERTIES CUDA_ARCHITECTURES 70)

add_subdirectory(test)