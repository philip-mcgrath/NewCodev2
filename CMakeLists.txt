cmake_minimum_required(VERSION 3.7)
project(NewCodev4)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp MultiPathProcessing.cpp transition.cpp functions.cpp random.cpp)
add_executable(NewCodev4 ${SOURCE_FILES})
find_package(GSL REQUIRED)
target_link_libraries(NewCodev4 GSL::gsl GSL::gslcblas)

include_directories(${GSL_INCLUDE_DIR})