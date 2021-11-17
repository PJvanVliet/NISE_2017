project(NISE_SAMPLES VERSION 1.0 LANGUAGES C)
add_subdirectory("tutorial")

add_custom_target(examples)
add_dependencies(examples Example::Tutorial)
