cmake_minimum_required(VERSION 3.14)
project(polyline-offset LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------- library -------
add_library(polyline_offset
    src/polyline_offset.cpp
)
target_include_directories(polyline_offset PUBLIC include)

# ------- tests ----------
enable_testing()

add_executable(test_offset test/test_offset.cpp)
target_link_libraries(test_offset PRIVATE polyline_offset)
add_test(NAME offset_tests COMMAND test_offset)
