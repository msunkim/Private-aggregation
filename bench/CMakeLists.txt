cmake_minimum_required(VERSION 3.16)
project(pnf_bench CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenFHE (OpenFHEConfig.cmake under /usr/local/lib/OpenFHE)
find_package(OpenFHE REQUIRED)

include_directories(${OpenFHE_INCLUDE})
include_directories(${OpenFHE_INCLUDE}/third-party/include)
include_directories(${OpenFHE_INCLUDE}/core)
include_directories(${OpenFHE_INCLUDE}/pke)
include_directories(${OpenFHE_INCLUDE}/binfhe)

link_directories(${OpenFHE_LIBDIR})

# Homebrew libomp (OpenFHE's flags point at a non-existent /usr/local/opt/llvm/lib)
link_directories(/usr/local/opt/libomp/lib)

# OpenFHE ships its flags as one space-separated string; apply them globally
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenFHE_CXX_FLAGS}")

# NTL / GMP for interpolation + root-finding over F_t
include_directories(/usr/local/include)

add_executable(threshold_spike threshold_spike.cpp)
target_link_libraries(threshold_spike ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(profile_b_bench profile_b_bench.cpp)
target_link_libraries(profile_b_bench ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

# Sigma-protocol prototype (NTL only; no OpenFHE dependency)
add_executable(sigma_spike sigma_spike.cpp)
target_link_libraries(sigma_spike ntl gmp)

add_executable(sigma_bench sigma_bench.cpp)
target_link_libraries(sigma_bench ntl gmp)

add_executable(ctsize ctsize.cpp)
target_link_libraries(ctsize ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(ksweep ksweep.cpp)
target_link_libraries(ksweep ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(probe_t probe_t.cpp)
target_link_libraries(probe_t ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(diag diag.cpp)
target_link_libraries(diag ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(thread_scale thread_scale.cpp)
target_link_libraries(thread_scale ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(profile_b128 profile_b128.cpp)
target_link_libraries(profile_b128 ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(ext_spike ext_spike.cpp)
target_link_libraries(ext_spike ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(ext_bench ext_bench.cpp)
target_link_libraries(ext_bench ${OpenFHE_SHARED_LIBRARIES} ntl gmp)

add_executable(ext_enc ext_enc.cpp)
target_link_libraries(ext_enc ${OpenFHE_SHARED_LIBRARIES} ntl gmp)
