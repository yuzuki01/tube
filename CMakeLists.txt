cmake_minimum_required(VERSION 3.25)

project(tube)

set(CMAKE_CXX_STANDARD 17)
set(PROJECT_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# OS type
if (WIN32)
    # if Windows
    message("OS is Windows")
    set(OS_INDICATOR "win")
    set(SHARED_LIBRARY_SUFFIX ".dll")
elseif (UNIX AND NOT APPLE)
    # if Linux
    message("OS is Linux")
    set(OS_INDICATOR "linux")
    set(SHARED_LIBRARY_SUFFIX ".so")
else ()
    message(FATAL_ERROR "unsupported OS.")
endif ()

set(LIBRARY_OUTPUT_PATH ${PROJECT_HOME_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_HOME_DIR}/bin)

include_directories(${PROJECT_HOME_DIR}/include)

file(GLOB SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
)

# 根据不同的目标执行不同的指令
add_executable(${PROJECT_NAME}
        ${SRC}
)

# check openMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP Found.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
