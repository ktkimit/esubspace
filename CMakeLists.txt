message(STATUS "Enriched Subspace Iteration Method")
message(STATUS "Cmake ${CMAKE_VERSION}\n")

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

#
# Project name, version, and main language
#
project(
  Esspace
  DESCRIPTION "Enriched subspace iteration method library"
  LANGUAGES Fortran
  )

file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION)
string(REPLACE "." ";" VERSION_LIST ${PROJECT_VERSION})
list(GET VERSION_LIST 0 PROJECT_VERSION_MAJOR)
list(GET VERSION_LIST 1 PROJECT_VERSION_MINOR)
list(GET VERSION_LIST 2 PROJECT_VERSION_PATCH)
unset(VERSION_LIST)

#
# Include packages
#
include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# include(CheckFortranSourceCompiles)
# include(CheckFortranSourceRuns)

#
# Set variables
#
set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")

#
# Create cmake files
#

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
  )

configure_package_config_file("${CMAKE_SOURCE_DIR}/cmake/config/Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  )

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  )

install(EXPORT "${TARGETS_EXPORT_NAME}"
  NAMESPACE "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  )

#
# Add subdirectories
#
add_subdirectory(src)

# if(BUILD_TESTING)
#   enable_testing()
# endif()
