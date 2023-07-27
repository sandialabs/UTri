# ############################################################################ #
#
#  UTri
#
# ############################################################################ #

# ---------------------------------------------------------------------------- #
# Set policy and CMake version requirements
# ---------------------------------------------------------------------------- #
cmake_minimum_required( VERSION 2.8 )


# ---------------------------------------------------------------------------- #
# Project definitions
# ---------------------------------------------------------------------------- #

# Enable Fortran to build the interface file
project(UTri C Fortran)

# Update the CMake module search path
set(CMAKE_MODULE_PATH ${UTri_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# Define the version
set(UTri_MAJOR_VERSION 0)
set(UTri_MINOR_VERSION 1)
set(UTri_PATCH_VERSION 1)
set(UTri_VERSION    ${UTri_MAJOR_VERSION}.${UTri_MINOR_VERSION})
if (UTri_PATCH_VERSION)
  set(UTri_VERSION "${UTri_VERSION}.${UTri_PATCH_VERSION}")
endif()

# ---------------------------------------------------------------------------- #
# External dependencies
# ---------------------------------------------------------------------------- #

# NetCDF : See FindNetCDF.cmake
find_package(NetCDF)

if ( NOT NETCDF_FOUND )
  message(FATAL_ERROR "Failed to locate NetCDF. Please define a NetCDF "
                      "installation with -D NETCDF_ROOT:PATH=<some path> or "
		      "define -D NETCDF_INCLUDE_DIR:PATH and "
		      "-D NETCDF_LIBRARY:FILEPATH explicitly")
endif()

# MPI
#  (1) Check if mpicc is in use. The check_include_file will
#      fail if not using a wrapper or mpi.h.
#  (2) If (1) fails, try to locate MPI_C_INCLUDE_PATH
#      either set explicitly with MPI_C_INCLUDE_PATH
#      or a call to find_package(MPI) if version > 2.8.6
include(CheckIncludeFile)
check_include_file(mpi.h MPI_H_IS_AVAILABLE)
if ( NOT MPI_H_IS_AVAILABLE )
  message(STATUS "Failed to include mpi.h")

  # The FindMPI.cmake before 2.8.6 was terrible
  # require MPI_C_INCLUDE_PATH to be set by the user
  if ( "${CMAKE_VERSION}" VERSION_LESS "2.8.6" )
    if ( NOT MPI_C_INCLUDE_PATH )
      message(WARNING "Using a cmake version (${CMAKE_VERSION}) with unreliable" 
                      " FindMPI.cmake. Please use newer version or define the MPI "
		      "include path explicitly with -D MPI_C_INCLUDE_PATH:PATH=")
    endif()
  else()
    find_package(MPI)
    if (NOT MPI_FOUND )
      message(WARNING "Could not locate MPI. Build will likely fail. "
	              "Please define MPI_C_INCLUDE_PATH explicitly.")
    endif()
  endif()

endif()  
        		    
# ---------------------------------------------------------------------------- #
# Fortran/C Interface
# ---------------------------------------------------------------------------- #
include(FortranCInterface)
FortranCInterface_HEADER(utri_FC.h
                         MACRO_NAMESPACE UTRIFC_
                         SYMBOLS utri_eos_handle_open
                                 utri_eos_handle_close
                                 utri_eos_status
                                 utri_eos_init
                                 utri_eos_delete
                                 utri_eos_isvinit
                                 utri_eos_isvquad
                                 utri_eos_pxdrv
                                 utri_eos_ptsat
                                 utri_eos_rtdrv
                                 utri_eos_redrv )

FortranCInterface_HEADER(utri_mpi_FC.h
                         MACRO_NAMESPACE UTRIMPIFC_
                         SYMBOLS utri_eos_handle_init_mpi_fortran )

if ( NOT FortranCInterface_GLOBAL_FOUND )
  message(WARNING "Could not determine the GLOBAL Fortran/C name mangling")
endif()

# ---------------------------------------------------------------------------- #
# Build the library
# ---------------------------------------------------------------------------- #

# Add the SOURCE_DIR and BINARY_DIR to the include paths
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

# Add NetCDF and MPI to the include paths
include_directories(${NETCDF_INCLUDE_DIRS})
include_directories(${MPI_C_INCLUDE_PATH})

# Library: utri (make utri)

# Source files
set(UTri_SOURCE_FILES
         eosdatapoints.c
	 utri_eos_mig.c
	 utri_eos_mpi.c
	 utri_eos_ncstatus.c
	 utri_eos_rptree.c
	 utri_eos_support.c)

# Header files       
set(UTri_HEADER_FILES
	 eosdatapoints.h
	 ncstatus.h
	 utri_eos_mig.h
	 utri_eos_mpi.h
	 utri_eos_support.h)

# Create the library target
add_library(utri ${UTri_SOURCE_FILES})

# Define the libraries needed to link against 
# this library
target_link_libraries(utri ${NETCDF_LIBRARIES} ${MPI_LIBRARIES})

# ---------------------------------------------------------------------------- #
# Installation 
# ---------------------------------------------------------------------------- #

# CMAKE_INSTALL_PREFIX defines the root install directory
# all paths are relative to this

# Include files
install(FILES ${UTri_HEADER_FILES} ${UTri_BINARY_DIR}/utri_FC.h
                                   ${UTri_BINARY_DIR}/utri_mpi_FC.h
        DESTINATION include)
     
# Library files
install(TARGETS utri
        EXPORT utri
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# Target file - must follow the INSTALL(TARGET utri ....) call
set(UTri_TARGET_FILE utri-targets.cmake)
install(EXPORT utri
        DESTINATION share/cmake/utri
        FILE ${UTri_TARGET_FILE})

# Config file, allow CMake projects to
#  find_package(UTri HINTS <install prefix>)
# and then UTri_INCLUDE_DIRS and UTri_LIBRARY are defined
set(UTri_CONFIG_FILE ${UTri_BINARY_DIR}/utri-config.cmake)
configure_file(${UTri_SOURCE_DIR}/utri-config.in 
               ${UTri_CONFIG_FILE} @ONLY)
install(FILES ${UTri_CONFIG_FILE}
        DESTINATION share/cmake/utri)
	     
# ---------------------------------------------------------------------------- #
# Distribution
# ---------------------------------------------------------------------------- #

# -- Binary package variables

# binary package name: default in the PROJECT name
set(CPACK_PACKAGE_NAME utri)

# Version information 
set(CPACK_PACKAGE_VERSION_MAJOR "${UTri_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${UTri_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${UTri_PATCH_LEVEL}")
set(CPACK_PACKAGE_VERSION       "${UTri_VERSION}")

# package file name
set(CPACK_PACKAGE_FILE_NAME utri-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME})

# Short description
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "UTri: Sandia Unstructured Triangle Tabular Interpolation Package")

# Vendor
set(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")



# --- Source package variables

# ZIP and GZIP'd and BZip2 TAR'd files, many other options
set(CPACK_SOURCE_GENERATOR "ZIP;TGZ;TBZ2")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "utri-${UTri_VERSION}")
  
# Files to ignore regular expressions are allowed
set(CPACK_SOURCE_IGNORE_FILES
     "^${UTri_SOURCE_DIR}/.hg"
     "^${UTri_SOURCE_DIR}/.cvs"
     "^${UTri_SOURCE_DIR}/.svn"
     ".*.swp"
     ".*[.]o"
     ".*~"
     "build"
     "CMakeFiles"
     "CMakeCache.txt"
     "cmake_install.cmake")

# This will add 'package' and 'package_source' targets
include(CPack)
