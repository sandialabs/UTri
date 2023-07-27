#include <stdio.h>
#include <stdlib.h>

#include "utri_eos_mpi.h"
#include "utri_eos_support.h"

/*
 * MPI broadcast function for the utri library
 *
 * inputs:
 */
void utri_eos_bcast_data_mpi( void * comm,
                              void * data,
                              int n,
                              int * error )
{
  *error = MPI_Bcast(data, n, MPI_CHAR, 0, *(MPI_Comm *)comm);
}

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
int utri_eos_handle_init_mpi( MPI_Comm comm )
{
  int rank;
  MPI_Comm * c;

  if (comm == MPI_COMM_NULL) {
    /* don't want to do MPI communications on a NULL communicator */
    return utri_eos_handle_init(NULL,NULL,1);
  }
  else {
    /* setup for MPI communications with rank 0 as master */
    rank = 0;
    if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
      /* call to MPI_Comm_rank failed */
      return -2;
    }

    /* use rank as master flag */
    if (rank == 0) rank = 1;
    else rank = 0;

    /* allocate space for an MPI_Comm */
    if ((c = (MPI_Comm *)malloc(sizeof(MPI_Comm))) == NULL) {
      perror("utri_eos_handle_mpi: unable to allocate memory for MPI_Comm");
      return -1;
    }

    /* final setup */
    *c = comm;
    return utri_eos_handle_init(utri_eos_bcast_data_mpi,c,rank);
  }
}

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
int utri_eos_handle_init_mpi_fortran( MPI_Fint * comm )
{
  return utri_eos_handle_init_mpi(MPI_Comm_f2c(*comm));
}
