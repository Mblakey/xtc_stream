#
# Michael Blakey 
# GROMACS XTC reverse engineering
# for streamed decompression
#
cmake_minimum_required(VERSION 3.5)
project(xtcstream C)

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

find_package(CURL REQUIRED)

include_directories(include)

set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set(SOURCES ${SOURCE_DIR}/xdrfile.c ${SOURCE_DIR}/xdrfile_trr.c ${SOURCE_DIR}/xdrfile_xtc.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3") 

add_executable(xtc-decompile ${SOURCE_DIR}/xtcdecompile.c ${SOURCES}) 
target_link_libraries(xtc-decompile CURL::libcurl)
