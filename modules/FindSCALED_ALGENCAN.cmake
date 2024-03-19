# - Try to find Scaled Algencan
# Once done this will define
#  SCALED_ALGENCAN_FOUND - System has Scaled Algencan
#  SCALED_ALGENCAN_LIB - The libraries needed to use Scaled Algencan

if (NOT SCALED_ALGENCAN_FOUND )

    find_path(SCALED_ALGENCAN_LIB
        NAMES libalgencan.a 
        PATHS "$ENV{SCALED_ALGENCAN_LIB}"
        "/opt/scaled-algencan/algencan-3.1.1_scaled/lib/"
        "~/scaled-algencan/algencan-3.1.1_scaled/lib/"
        )

    include(FindPackageHandleStandardArgs)

    # handle the QUIETLY and REQUIRED arguments and set SCALED_ALGENCAN_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(SCALED_ALGENCAN DEFAULT_MSG SCALED_ALGENCAN_LIB)

    mark_as_advanced(SCALED_ALGENCAN_LIB)

endif(NOT SCALED_ALGENCAN_FOUND)

