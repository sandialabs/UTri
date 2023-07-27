/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of UTri, substitue the full text  */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"

#include "utri_eos_support.h"
#include "ncstatus.h"

/*
 * check the netcdf status and if needed print out an 80 character
 * wide, multi-line error message using the MIG printing utilities
 *
 */
int bad_netcdf_status(int status,char *ncfun,const char *database,char *fun,char *file,int line)
{
  int length;
  char *s;
  char m[81];
  int tot;

  if (status != 0) {
    /* overestimate length */
    length = 13+15+strlen(nc_strerror(status))+10+strlen(ncfun)+22+strlen(database)+12+strlen(fun)+8+strlen(file)+8+15+2;
    if ((s = (char *)malloc(sizeof(char)*length)) == NULL) {
      reportUTriError("bad_netcdf_status: failed to allocate memory");
      return status;
    }
    tot = snprintf(s,length,
            "NetCDF error %d \"%s\" calling: %s, database: %s, function: %s, file: %s, line: %d",
	     status,nc_strerror(status),ncfun,database,fun,file,line);

    if (tot < 1) {
      reportUTriError("bad_netcdf_status: failed to generate NetCDF error string");
      return status;
    }
    if (tot >= length) tot = length-1;
    length = 0;
    while (length < tot) {
      strncpy(m,&s[length],80);
      m[80] = '\0';
      reportUTriError(m);
      length += 80;
    }
    free(s);
  }
  
  return status;
}
