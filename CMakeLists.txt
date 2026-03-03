cmake_minimum_required(VERSION 3.14)
project(polyline_offset VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library
add_library(polyline_offset
    src/polyline_offset.cpp
)
target_include_directories(polyline_offset PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Tests
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_executable(test_offset tests/test_offset.cpp)
    target_link_libraries(test_offset PRIVATE polyline_offset)
    add_test(NAME test_offset COMMAND test_offset)
endif()
