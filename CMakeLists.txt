cmake_minimum_required(VERSION 3.1)
project(Bezier
  LANGUAGES CXX
  VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-fPIC)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

set(Bezier_SRC
  src/bezier.cpp
  src/polycurve.cpp
  )

set(Bezier_INC
  include/Bezier/declarations.h
  include/Bezier/legendre_gauss.h
  include/Bezier/bezier.h
  include/Bezier/polycurve.h
  )

# Options

if(NOT LEGENDRE_GAUSS_PRECISION)
  message(STATUS "Setting Legendre Gauss precision to default value (10).")
  set(LEGENDRE_GAUSS_PRECISION  10)
endif()

option(BUILD_SHARED_LIBS "Build shared library (.so) instead of static one (.a)" OFF)

# build rules

if(BUILD_SHARED_LIBS)
  add_library(bezier SHARED ${Bezier_SRC})
else()
  add_library(bezier STATIC ${Bezier_SRC})
endif()

target_include_directories(bezier PUBLIC
  $<BUILD_INTERFACE:${Bezier_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_compile_definitions(bezier PRIVATE LEGENDRE_GAUSS_N=${LEGENDRE_GAUSS_PRECISION})

set_target_properties(bezier PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(bezier PROPERTIES PUBLIC_HEADER "${Bezier_INC}")

# install rules
install(TARGETS bezier
  EXPORT bezier-export DESTINATION "lib"
  PUBLIC_HEADER DESTINATION "include/Bezier")
install(EXPORT bezier-export DESTINATION "lib/cmake/Bezier" FILE BezierConfig.cmake)
