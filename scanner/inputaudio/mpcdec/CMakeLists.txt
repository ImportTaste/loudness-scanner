cmake_minimum_required(VERSION 2.8.9)

include(utils)

find_package(MPCDEC)
find_package(PkgConfig)
pkg_check_modules(GMODULE20 gmodule-2.0)

set(SUMMARY_MPCDEC_FOUND ${MPCDEC_FOUND} CACHE INTERNAL "")

if(GMODULE20_FOUND AND MPCDEC_FOUND AND
   INPUT_INCLUDE_DIR AND EBUR128_INCLUDE_DIR AND
   NOT DISABLE_GLIB20 AND NOT DISABLE_MPCDEC)
  include_directories(SYSTEM ${MPCDEC_INCLUDE_DIRS})
  include_directories(${INPUT_INCLUDE_DIR} ${EBUR128_INCLUDE_DIR})

  add_library(input_musepack MODULE input_mpcdec.c)

  target_link_libraries(input_musepack ${MPCDEC_LIBRARIES})
  install(TARGETS input_musepack DESTINATION ${CMAKE_INSTALL_LIBDIR})

  list(APPEND INPUT_MPCDEC_CFLAGS  ${GMODULE20_CFLAGS})
  list(APPEND INPUT_MPCDEC_LDFLAGS ${GMODULE20_LDFLAGS})

  to_space_list(INPUT_MPCDEC_CFLAGS)
  to_space_list(INPUT_MPCDEC_LDFLAGS)

  set_target_properties(input_musepack PROPERTIES
          COMPILE_FLAGS ${INPUT_MPCDEC_CFLAGS}
             LINK_FLAGS ${INPUT_MPCDEC_LDFLAGS})
endif()
