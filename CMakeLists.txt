cmake_minimum_required(VERSION 3.7)
project(NewCodev2)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp MultiPathProcessing.cpp transition.cpp functions.cpp)
add_executable(NewCodev2 ${SOURCE_FILES})