cmake_minimum_required (VERSION 3.14)
project (FrameSim)

set (SOURCES circuit.cpp simulator.cpp noise.cpp nonsparsesim.cpp simulator_base.cpp)

add_definitions(-DNO_THREADS)
add_definitions(-DCMAKE_CXX_FLAGS="-Werror -Wall -Wextra")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (LINUX)
        add_definitions(-D_GLIBCXX_DEBUG)
    endif()
endif()

add_library(FrameSim SHARED ${SOURCES})
