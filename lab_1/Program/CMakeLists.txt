cmake_minimum_required(VERSION 3.17)
project(Spline)

set(CMAKE_CXX_STANDARD 14)

add_executable(Spline First/main.cpp First/Function.cpp First/Function.h)
add_executable(lagr First/main.cpp)
