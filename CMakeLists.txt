cmake_minimum_required(VERSION 3.16)
project(ScaledAlgencan-CMake)

##### Include Dirs 

include_directories(
    SYSTEM
    includes
    "${PROJECT_BINARY_DIR}"
    )

##### Modules folder

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/modules)

##### Scaped Algencan setup

find_package (SCALED_ALGENCAN)

add_library(algencan STATIC IMPORTED)
set_target_properties(algencan PROPERTIES
    IMPORTED_LOCATION "${SCALED_ALGENCAN_LIB}/libalgencan.a"
    INTERFACE_LINK_LIBRARIES "gfortran;m"
    )
set_target_properties(algencan PROPERTIES
    IMPORTED_LOCATION "${SCALED_ALGENCAN_LIB}/libalgencan.so"
    )


##### Compile options 

add_compile_options(
    -std=c++20
    # ALGENCAN
    -lalgencan 
    -lgfortran 
    -lm 
    -lopenblas 
    -zmuldefs 
    )

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(main ${SOURCES})

target_link_libraries(main PRIVATE algencan)
