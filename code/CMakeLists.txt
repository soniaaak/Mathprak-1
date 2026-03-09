cmake_minimum_required(VERSION 3.15)
project(first_result)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(result result.cpp)
add_executable(result_simple result_simple.cpp)
add_executable(gear gear_mesh.cpp)

find_path(GMSH_INCLUDE_DIR gmsh.h)
find_library(GMSH_LIBRARY gmsh)

if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIBRARY)
    message(FATAL_ERROR "gmsh not found. Install it with: brew install gmsh")
endif()

target_include_directories(result PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(result_simple PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(gear PRIVATE ${GMSH_INCLUDE_DIR})
target_link_libraries(result PRIVATE ${GMSH_LIBRARY})
target_link_libraries(result_simple PRIVATE ${GMSH_LIBRARY})
target_link_libraries(gear PRIVATE ${GMSH_LIBRARY})