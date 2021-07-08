cmake_minimum_required(VERSION 3.17)
project(Pendulum)

set(CMAKE_CXX_STANDARD 14)

add_executable(Pendulum main.cpp quanc8.cpp quanc8.h rkf45.cpp rkf45.h)
