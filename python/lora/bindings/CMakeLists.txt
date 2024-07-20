# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

########################################################################
# Check if there is C++ code at all
########################################################################
if(NOT lora_sources)
    message(STATUS "No C++ sources... skipping python bindings")
    return()
endif(NOT lora_sources)

########################################################################
# Check for pygccxml
########################################################################
gr_python_check_module_raw("pygccxml" "import pygccxml" PYGCCXML_FOUND)

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND lora_python_files
        decode_python.cc
        demod_python.cc
        encode_python.cc
        mod_python.cc
        pyramid_demod_python.cc
        weak_demod_python.cc
        python_bindings.cc)

gr_pybind_make_oot(lora ../../.. gr::lora "${lora_python_files}")

# copy bindings extension for use in QA test module
add_custom_command(
    TARGET lora_python
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:lora_python>
            ${PROJECT_BINARY_DIR}/test_modules/gnuradio/lora/)

install(
    TARGETS lora_python
    DESTINATION ${GR_PYTHON_DIR}/gnuradio/lora
    COMPONENT pythonapi)
