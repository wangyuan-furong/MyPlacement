# cmake needs this line
cmake_minimum_required(VERSION 3.5)

# Define project name
project(myplacement)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CImg 编译时禁用显示（本项目仅用文件 IO 画图；系统无 X11 头文件，不能链接 -lX11）
add_definitions(-Dcimg_display=0)

# Declare the executable target built from your sources
add_executable(
    myplace
    myplace_main.cpp
    PlaceData/placedata.cpp
    MyPlace/myplace.cpp
    Plot/plot.cpp
)

include_directories(
    ${PROJECT_SOURCE_DIR}/MyPlace
    ${PROJECT_SOURCE_DIR}/PlaceData
    ${PROJECT_SOURCE_DIR}/Plot
    ${PROJECT_SOURCE_DIR}/Common
    ${PROJECT_SOURCE_DIR}/Library/Eigen
    ${PROJECT_SOURCE_DIR}/Library/CImg
)

target_link_libraries(myplace "-lpthread -lm -ldl")
