cmake_minimum_required(VERSION 3.0)
project(nearlab_utils)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11 -g -Wall)

find_package(catkin REQUIRED COMPONENTS
  roscpp
)

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
set(EIGEN3_INCLUDE_DIRS /usr/share/eigen)

#########
# BUILD #
#########
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIRS}
)


add_library(nearlab_utils SHARED
               src/orbitPropagator.cpp
               src/quatMath.cpp
               src/rungeKutta.cpp
	       src/attitudePropagator.cpp)

target_link_libraries(nearlab_utils ${catkin_LIBRARIES})

#set_target_properties(nearlab_utils PROPERTIES PUBLIC_HEADER include/*.h)


catkin_package(
  CATKIN_DEPENDS roscpp
  INCLUDE_DIRS include
  LIBRARIES nearlab_utils
)

install(DIRECTORY include/${PROJECT_NAME}/
	DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})

install( TARGETS nearlab_utils
         DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}/${PROJECT_NAME} )
