add_executable(${ExampleName} ${ExampleName}.F90)

target_link_libraries(${ExampleName} PRIVATE opencmisslibs)

install(TARGETS ${ExampleName} DESTINATION bin)

configure_file(visualise.cmgui visualise.cmgui)

configure_file(inputs/lvq-842.in.gpactiv ${PROJECT_BINARY_DIR}/lvq-842.in.gpactiv)
configure_file(inputs/lvq-842.in ${PROJECT_BINARY_DIR}/lvq-842.in)