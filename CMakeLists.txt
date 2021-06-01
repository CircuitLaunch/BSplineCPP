cmake_minimum_required(VERSION 3.10)

project(BSplineCPP VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)

set(PROGRAMS test)
set(CORELIBS m)

configure_file(BSplineCPPConfig.h.in BSplineCPPConfig.h)

foreach(program ${PROGRAMS})
  add_executable(${program} test.cpp)
  target_include_directories(${program} PUBLIC "${PROJECT_BINARY_DIR}")
  target_link_libraries(${program} ${CORELIBS} BSpline)
endforeach(program)
