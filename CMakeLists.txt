cmake_minimum_required(VERSION 4.0)
project(MAIN)

# 添加这行：强制使用 C99 标准
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(inc)
set(SRC_LIST src/main.c src/algebra.c)
add_executable(main ${SRC_LIST})