cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(Example_frame2d
  LANGUAGES Fortran)

find_package(Esspace REQUIRED)

add_executable(${PROJECT_NAME} esspace_main.f90)

target_link_libraries(${PROJECT_NAME} PRIVATE Esspace::Esspace)

file(GLOB INPUT_DATA "frame2d_ex.*")
file(COPY ${INPUT_DATA}
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
