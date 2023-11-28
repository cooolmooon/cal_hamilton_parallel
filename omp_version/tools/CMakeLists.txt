add_library(input STATIC input.cpp)
add_library(timer STATIC timer.cpp)
add_library(spline STATIC spline.cpp)
add_library(diago STATIC diago.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(timer PUBLIC OpenMP::OpenMP_CXX)
endif()

find_package(LAPACK REQUIRED)

target_link_libraries(diago lapack lapacke)