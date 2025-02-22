cmake_minimum_required (VERSION 2.8)
project (edgepart)

set(CMAKE_VERBOSE_MAKEFILE off)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package (Gflags REQUIRED)
find_package (Glog REQUIRED)
find_package (Boost REQUIRED)

# # SET a default build type if none was specified
# if(NOT CMAKE_BUILD_TYPE)
#   message(STATUS "Setting build type to 'Release' as none was specified.")
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
#   # SET the possible values of build type for cmake-gui
# endif()
# set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -pthread -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_EXE_LINK_FLAGS_RELEASE} -Ofast -DNDEBUG")

# build a (static) library for the exposed API
add_library(edgepart src/edgepart.cpp src/nv_partitioner.cpp)
target_include_directories(edgepart PUBLIC include)
target_link_libraries (edgepart glog)

if (NOT EDGEPART_LIB_ONLY)
  add_subdirectory(threadpool11)
  add_executable (main
      src/main.cpp
      src/util.cpp
      src/sort.cpp
      src/graph.cpp
      src/nv_partitioner.cpp
      src/sne_partitioner.cpp
      src/random_partitioner.cpp
      src/hsfc_partitioner.cpp
      src/dbh_partitioner.cpp
      src/conversions.cpp
      src/shuffler.cpp)
  add_executable (graph2edgelist
      src/graph2edgelist.cpp
      src/util.cpp
      src/conversions.cpp)

  target_include_directories(main PRIVATE include ./threadpool11/include)
  target_include_directories(graph2edgelist PRIVATE include)

  target_link_libraries (main glog gflags threadpool11 edgepart)
  target_link_libraries (graph2edgelist glog gflags)
endif()
