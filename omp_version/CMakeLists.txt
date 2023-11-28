cmake_minimum_required(VERSION 3.9)
project(big)
add_subdirectory(tools)

add_executable(${PROJECT_NAME} main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(big PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE input)
target_link_libraries(${PROJECT_NAME} PRIVATE timer)
target_link_libraries(${PROJECT_NAME} PRIVATE spline)
target_link_libraries(${PROJECT_NAME} PRIVATE diago)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/tools)