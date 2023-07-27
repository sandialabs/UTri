/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of UTri, substitue the full text  */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef UTRI_EOS_MPI_H
#define UTRI_EOS_MPI_H

#include "utri_mpi_FC.h"
#include "mpi.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
 * Note that there are two init functions in the UTri EOS MPI
 * interface. The first for C/C++ and the second for FORTRAN.
 *
 */

/*
 * Initialize the utri eos library with the given MPI communicator.
 *
 * inputs:
 *
 *   comm -- a valid MPI communicator, can be MPI_COMM_NULL
 *
 * return:
 *
 *   integer handle for this instance of the utri eos library
 */
int utri_eos_handle_init_mpi( MPI_Comm comm );

/*
 * Initialize the utri eos library with the given MPI
 * communicator. This FORTRAN version converts the fortran
 * communicator to a c version. FORTRAN linking is assumed, so
 * variables are passed by reference.
 *
 * inputs:
 *
 *   comm -- a valid FORTRAN MPI communicator, can be MPI_COMM_NULL equivalent
 *
 * return:
 *
 *   integer handle for this instance of the utri eos library
 */
int utri_eos_handle_init_mpi_fortran( MPI_Fint * comm );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
