cmake_minimum_required(VERSION 3.15)
project(laba2_do)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(gear_anim result2.cpp)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
target_link_libraries(gear_anim PRIVATE ${VTK_LIBRARIES})

find_path(GMSH_INCLUDE_DIR gmsh.h)
find_library(GMSH_LIBRARY gmsh)

if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIBRARY)
    message(FATAL_ERROR "gmsh not found. Install it with: brew install gmsh")
endif()

target_include_directories(gear_anim PRIVATE ${GMSH_INCLUDE_DIR})
target_link_libraries(gear_anim PRIVATE ${GMSH_LIBRARY})
