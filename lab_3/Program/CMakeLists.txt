cmake_minimum_required(VERSION 3.17)
project(Differencial)

set(CMAKE_CXX_STANDARD 14)

add_executable(Differencial1 main.cpp rkf45.cpp rkf45.h)