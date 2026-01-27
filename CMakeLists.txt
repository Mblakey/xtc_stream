cmake_minimum_required(VERSION 3.5)
project(xtcstream C)

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

include(CheckCCompilerFlag)
check_c_compiler_flag("-msse2" COMPILER_SUPPORTS_SSE2)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3") 

include_directories(include)
find_package(CURL REQUIRED)

set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(ORIG_SOURCES ${SOURCE_DIR}/xdrfile.c ${SOURCE_DIR}/xdrfile_trr.c ${SOURCE_DIR}/xdrfile_xtc.c)

set(OPT_SOURCES ${SOURCE_DIR}/xdrfile_opt.c ${SOURCE_DIR}/xdrfile_trr.c ${SOURCE_DIR}/xdrfile_xtc.c)

add_executable(xtc-decompile ${SOURCE_DIR}/xtcdecompile.c ${ORIG_SOURCES}) 
add_executable(xtc-decompile-opt ${SOURCE_DIR}/xtcdecompile.c ${OPT_SOURCES}) 
add_executable(test-harness ${SOURCE_DIR}/test_harness.c) 
add_executable(mkbitjmptlb ${SOURCE_DIR}/mkjmptlb.c) 

if (COMPILER_SUPPORTS_SSE2)
  target_compile_options(xtc-decompile PRIVATE -msse2)
  target_compile_options(xtc-decompile-opt PRIVATE -msse2)
  target_compile_options(test-harness PRIVATE -msse2)
endif()

target_link_libraries(xtc-decompile CURL::libcurl)
target_link_libraries(xtc-decompile-opt CURL::libcurl)
