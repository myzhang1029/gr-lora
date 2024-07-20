find_package(PkgConfig)

PKG_CHECK_MODULES(PC_GR_LORA gnuradio-lora)

FIND_PATH(
    GR_LORA_INCLUDE_DIRS
    NAMES gnuradio/lora/api.h
    HINTS $ENV{LORA_DIR}/include
        ${PC_LORA_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    GR_LORA_LIBRARIES
    NAMES gnuradio-lora
    HINTS $ENV{LORA_DIR}/lib
        ${PC_LORA_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-loraTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GR_LORA DEFAULT_MSG GR_LORA_LIBRARIES GR_LORA_INCLUDE_DIRS)
MARK_AS_ADVANCED(GR_LORA_LIBRARIES GR_LORA_INCLUDE_DIRS)
