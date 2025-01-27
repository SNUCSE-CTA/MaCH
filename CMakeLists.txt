cmake_minimum_required(VERSION 3.13)
project(GraphLib)

set(CMAKE_CXX_STANDARD 20)
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")

include_directories(lib)

add_executable(MaCH driver/subhypergraph-matching.cc)
add_executable(HypergraphPreprocessing driver/hypergraph-preprocessing.cc)
