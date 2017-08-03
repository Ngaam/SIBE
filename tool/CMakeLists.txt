# Collect source files
file(GLOB_RECURSE srcs ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Build each source file independently
foreach(source ${srcs})
  get_filename_component(name ${source} NAME_WE)

  # sibe target already exits
  if(name MATCHES "sibe")
    set(name ${name}.bin)
  endif()

  # target
  add_executable(${name} ${source})
  target_link_libraries(${name} ${Sibe_LINK})
  sibe_default_properties(${name})

  # set back RUNTIME_OUTPUT_DIRECTORY
  sibe_set_runtime_directory(${name} "${PROJECT_BINARY_DIR}/tool")
  sibe_set_solution_folder(${name} tool)

  # restore output name without suffix
  if(name MATCHES "sibe.bin")
    set_target_properties(${name} PROPERTIES OUTPUT_NAME sibe)
  endif()

  # Install
  install(TARGETS ${name} DESTINATION bin)
endforeach(source)
