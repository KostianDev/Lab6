cmake_minimum_required(VERSION 3.26)
project(Lab6 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)

set(SDL2_PATH "C:/Users/user/OneDrive/Documents/SDL2-2.30.2/x86_64-w64-mingw32")
set(SDL2TTF "C:/Users/user/OneDrive/Documents/SDL2_ttf-2.22.0/x86_64-w64-mingw32")

find_package(SDL2 REQUIRED)
find_package(SDL2TTF REQUIRED)
include_directories(${SDL2_INCLUDE_DIR} ${SDL2TTF_INCLUDE_DIR})

add_executable(Lab6 Lab6.c)

target_link_libraries(${PROJECT_NAME} ${SDL2_LIBRARY} ${SDL2TTF_LIBRARY})