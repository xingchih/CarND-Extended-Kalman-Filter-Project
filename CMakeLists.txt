project(ExtendedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp src/tools.cpp src/FusionEKF.cpp src/kalman_filter.cpp)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
#include_directories(/usr/local/opt/openssl/include)
include_directories(/usr/local/Cellar/openssl/1.0.2k/include)
include_directories(/Users/xingchihe/Documents/udacity/CarND/Term2/CarND-Extended-Kalman-Filter-Project/uWebSockets/src)
link_directories(/usr/local/lib)
#link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/openssl/1.0.2k/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)
link_directories(/Users/xingchihe/Documents/udacity/CarND/Term2/CarND-Extended-Kalman-Filter-Project/uWebSockets/build)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 


add_executable(ExtendedKF ${sources})

target_link_libraries(ExtendedKF z ssl uv uWS)
