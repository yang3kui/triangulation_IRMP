cmake_minimum_required(VERSION 2.8.3)
project(triangulation)

set(CMAKE_BUILD_TYPE RELEASE)
add_definitions(-std=c++11)
add_compile_options(-O3)

find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(
	${CERES_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
	)

file(GLOB_RECURSE CPP_FILES src/*.cpp)
file(GLOB_RECURSE CPP_HEADERS src/*.h)


add_executable(triangulation ${CPP_FILES} ${CPP_HEADERS})

target_link_libraries(triangulation ${CERES_LIBRARIES})
