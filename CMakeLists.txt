cmake_minimum_required(VERSION 2.8.7)
if(POLICY CMP0046)
  cmake_policy(SET CMP0046 NEW)
endif()
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

# ---[ Sibe project
project(Sibe C CXX)

# ---[ Sibe version
set(SIBE_TARGET_VERSION "1.3.1" CACHE STRING "Sibe logical version")
set(SIBE_TARGET_SOVERSION "1.4.1" CACHE STRING "Sibe soname version")
add_definitions(-DSIBE_VERSION=${SIBE_TARGET_VERSION})


# for errors in /usr/include/string.h
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

if(CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
else ()
  set(CMAKE_CXX_STANDARD 11)
endif ()


# ---[ Using cmake scripts and modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

#include(ExternalProject)

include(cmake/Utils.cmake)
include(cmake/Targets.cmake)
include(cmake/Misc.cmake)
include(cmake/Summary.cmake)
include(cmake/ConfigGen.cmake)

# ---[ Options
sibe_option(CPU_ONLY  "Build Sibe without CUDA support" OFF)
#sibe_option(USE_CUDNN "Build Sibe with cuDNN library support" ON IF NOT CPU_ONLY)
sibe_option(BUILD_SHARED_LIBS "Build shared libraries" ON)
sibe_option(BUILD_docs   "Build documentation" OFF IF UNIX OR APPLE)
sibe_option(USE_LEVELDB "Build with levelDB" OFF)
sibe_option(USE_LMDB "Build with lmdb" OFF)

# ---[ Dependencies
include(cmake/Dependencies.cmake)

# ---[ Flags
if(UNIX OR APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall")
endif()

sibe_set_sibe_link()

#if(USE_libstdcpp)
##  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ -std=c++11")
#  message("-- Warning: forcing libstdc++ (controlled by USE_libstdcpp option in cmake)")
#endif()

if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
else ()
  set(CMAKE_CXX_STANDARD 11)
endif ()

add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)

# ---[ Warnings
sibe_warnings_disable(CMAKE_CXX_FLAGS -Wno-sign-compare -Wno-uninitialized)

# ---[ Config generation
configure_file(cmake/Templates/sibe_config.h.in "${PROJECT_BINARY_DIR}/sibe_config.h")

# ---[ Includes
set(Sibe_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories(${Sibe_INCLUDE_DIR} ${PROJECT_BINARY_DIR})
include_directories(BEFORE src) # This is needed for gtest.

# ---[ Subdirectories
#add_subdirectory(src/gtest)
add_subdirectory(src/sibe)
add_subdirectory(tool)
#add_subdirectory(examples)
#add_subdirectory(docs)

# ---[ Linter target
add_custom_target(lint COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/cmake/lint.cmake)

# ---[ Configuration summary
sibe_print_configuration_summary()

# ---[ Export configs generation
sibe_generate_export_configs()

sibe_extract_sibe_version()
