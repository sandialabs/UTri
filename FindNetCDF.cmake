include(FindPackageHandleStandardArgs)

if ( NOT NETCDF_INCLUDE_DIRS )

  find_path(NETCDF_INCLUDE_DIR netcdf.h 
            HINTS ${NETCDF_ROOT} ENV NETCDF_ROOT
	    PATH_SUFFIXES include)

  if ( NETCDF_INCLUDE_DIR )
    set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
  endif()

endif()

if ( NOT NETCDF_LIBRARIES )

  find_library(NETCDF_LIBRARY netcdf 
               HINTS ${NETCDF_ROOT} ENV NETCDF_ROOT
	       PATH_SUFFIXES lib)

  if ( NETCDF_LIBRARY )	     
    set(NETCDF_LIBRARIES ${NETCDF_LIBRARY})
  endif()

endif()

find_package_handle_standard_args(NetCDF DEFAULT_MSG
				  NETCDF_LIBRARIES
                                  NETCDF_INCLUDE_DIRS)


            
