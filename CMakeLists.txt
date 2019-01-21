cmake_minimum_required(VERSION 3.0)
project(nearlab_utils)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11 -g -Wall)

find_package(catkin REQUIRED COMPONENTS
  roscpp
)

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
set(EIGEN3_INCLUDE_DIRS /usr/share/eigen)

catkin_package(CATKIN_DEPENDS
  roscpp
)

#########
# BUILD #
#########
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIRS}
)


add_library(nearlab_utils
               src/orbitPropagator.cpp
               src/quatMath.cpp
               src/rungeKutta.cpp)

set_target_properties(nearlab_utils PROPERTIES PUBLIC_HEADER include/*.h)


# install(TARGETS mylib
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
