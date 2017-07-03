cmake_minimum_required(VERSION 3.7)
project(untitled1)

set(CMAKE_C_STANDARD 99)

include_directories(src/AccelAndMagnet src/libfixmath src/libfixmatrix src/libcsv)
link_directories(src/AccelAndMagnet src/libfixmath src/libfixmatrix src/libcsv)

file(GLOB_RECURSE SOURCE_FILES ./src/*.c)

add_executable(untitled1 ${SOURCE_FILES})