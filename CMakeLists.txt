cmake_minimum_required(VERSION 3.27)
project(rscoder)

set(CMAKE_CXX_STANDARD 23)

add_library(rscoder INTERFACE)

target_include_directories(rscoder INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)

find_package(TBB REQUIRED)
target_link_libraries(rscoder INTERFACE TBB::tbb)
target_compile_features(rscoder INTERFACE cxx_std_23)