cmake_minimum_required(VERSION 3.20)
project(CelestialSPHSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(celestial_core
    src/physics/gravity/BarnesHutOctree.cpp
)

target_include_directories(celestial_core PUBLIC include)

target_compile_options(celestial_core PRIVATE -O3 -march=native)

add_executable(octree_demo src/app/main.cpp)
target_link_libraries(octree_demo PRIVATE celestial_core)
