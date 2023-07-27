/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of UTri, substitue the full text  */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

/* ncstatus.h */

#ifndef NC_STATUS_H
#define NC_STATUS_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
   This function is called by the netcdf status check macro to print
   an error message if status is non-zero. This may be implemented
   differently depending on how the output is desired. The status is
   returned if control does not otherwise terminate inside this
   function.
*/

int bad_netcdf_status( int status,
		       char *ncfun,
		       const char *database,
		       char *fun,
		       char *file,
		       int line );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
