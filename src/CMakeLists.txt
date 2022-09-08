set(src_dir 
  linalg 
  )

foreach(directory ${src_dir})
  add_subdirectory(${directory})
endforeach()

set(src
  esspace.f90
  )

add_library(${PROJECT_NAME}
  ${src}
  )

target_link_libraries(${PROJECT_NAME}
  PUBLIC
  Linalg
  # "$<BUILD_INTERFACE:Linalg>"
  )

# A module directory is generated and mod files are put there.
set(MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod)
set_target_properties(${PROJECT_NAME}
  PROPERTIES
  Fortran_MODULE_DIRECTORY
  ${MOD_DIR}
  )

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${MOD_DIR}>
  $<INSTALL_INTERFACE:mod>
  )

install(TARGETS ${PROJECT_NAME}
  EXPORT ${TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

install(DIRECTORY 
  ${MOD_DIR} 
  DESTINATION 
  "${CMAKE_INSTALL_PREFIX}"
  )
