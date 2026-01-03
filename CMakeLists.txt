cmake_minimum_required(VERSION 3.5)
project(xtcstream C)

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

include(CheckCCompilerFlag)
check_c_compiler_flag("-msse2" COMPILER_SUPPORTS_SSE2)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3") 

include_directories(include)
find_package(CURL REQUIRED)

set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(SOURCES ${SOURCE_DIR}/xdrfile.c ${SOURCE_DIR}/xdrfile_trr.c ${SOURCE_DIR}/xdrfile_xtc.c)

add_executable(xtc-decompile ${SOURCE_DIR}/xtcdecompile.c ${SOURCES}) 

if (COMPILER_SUPPORTS_SSE2)
  target_compile_options(xtc-decompile PRIVATE -msse2)
endif()

target_link_libraries(xtc-decompile CURL::libcurl)
