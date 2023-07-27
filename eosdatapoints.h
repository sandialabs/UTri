/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of UTri, substitue the full text  */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

/* eosdatapoints.h*/

#ifndef EOSDATAPOINTS_H
#define EOSDATAPOINTS_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
   The eos_data_points structure is used to save tabular EOS data in
   array fashion to support vectorization of interpolation routines
  
   The 'n' entry gives the number of points present in each of the
   arrays. The array data type is commented below
*/
typedef struct eos_data_points
{
  int n;         /* number of points in each of the following arrays */
  int *phase;    /* phase of material at this state */
  double *rho;   /* density */
  double *T;     /* temperature */
  double *p;     /* pressure */
  double *E;     /* specific internal energy */
  double *S;     /* specific entropy */
  double *F;     /* specific Helmholtz free energy */
  double *G;     /* specific Gibb's free energy */
  double *H;     /* specific enthalpy */
  double *dpdr;  /* rho derivative of p at constant T */
  double *dpdT;  /* T derivative of p at constant rho */
  double *dEdr;  /* rho derivative of E at constant T */
  double *dEdT;  /* T derivative of E at constant rho */
  double *cs;    /* adiabatic sound speed */
  double *dcsdr; /* rho derivative of cs at constant S */
  double *visc;  /* viscosity */
  double *thcon; /* thermal conductivity */
} eos_data_points_t;

/*
   Method used to (re)allocate the arrays in the eos_data_points_t
   structure to the size 'l'. Calling/return semantics are similar to
   the standard C library realloc() function. One should not pass l=0
   to this function to attempt to free the resources associated with
   an eos_data_points_t structure. Use the eospointsfree function
   instead.
*/
void *eospointsrealloc(void *p,int l);

/*
   Method to free resources associated with the eos_data_points_t
   structure. The pointer p must have been initialized by a call to
   eospointsrealloc.
*/
void eospointsfree(void *p);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
