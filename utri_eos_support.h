/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of UTri, substitue the full text  */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

/* utri_eos_support.h*/

#ifndef UTRI_EOS_SUPPORT_H
#define UTRI_EOS_SUPPORT_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include "eosdatapoints.h"

/*
 * Error reporting function for utri interpolation code. The function
 * appends the message to the current error string.
 *
 * inputs:
 *
 *   message -- null terminated character string
 */
void reportUTriError( const char * message );

/*
 * return the error string and clear the internal buffer
 *
 * Caller must free memory.
 */
char * utri_error_string_get();

/*
 * Initialize a new handle.
 *
 * inputs:
 *
 *   fun    -- pointer to broadcast function
 *   data   -- opaque pointer for broadcast function data
 *   master -- flag for controlling process, non-zero means master
 *
 * return:
 *
 *   integer handle or -1 on error.
 */
int utri_eos_handle_init( void (*fun)( void *, void *, int, int * ),
                          void * data,
                          int master );

/*
 * Close and free resources associated with a utri eos handle.
 *
 * inputs:
 *
 *   handle -- integer handle to utri data
 */
void utri_eos_handle_free( const int handle );

/*
   Structure to represent a triangle in an unstructured grid and
   support evaluation of the barycentric coordinates. Triangle
   vertices are given in anti-clockwise orientation.
*/
typedef struct eos_table_tri
{
  double a[6]; /* precomputed constants for barycentric coordinate calculation*/
  int p[3];    /* triangle vertex index to the data points array*/
} eos_table_tri_t;

/*
   Structure to store a RP tree for fast lookup of triangles in an
   unstructured grid on a rectangular domain.
  
   The tree is generated through a series of splits alternating first
   on the x and then y directions. There are 2^slevels leaf nodes at
   the bottom of the tree. Each of these leaf nodes has a pair of
   integers stored consecutively in the trilookup array. The first
   integer is an index into the trilist array. The second is the
   number of triangles in the leaf. The slice of the trilist array
   defined by this pair contains the triangle numbers for the
   triangles in the leaf. These numbers index into the triangle array
   associated with the tree.
*/
typedef struct eos_table_rptree
{
  double xbounds[2]; /* x values for overall bounding rectangle */
  double ybounds[2]; /* y values for overall bounding rectangle */
  double *splits;    /* array of x,y split values */
  int *trilookup;    /* array to lookup triangles in a split */
  int *trilist;      /* array of triangles in splits */
  int slevels;       /* number of splits */
} eos_table_rptree_t;

/*
   Structure stored in memory and on disk for holding an unstructured
   triangular EOS table with PCA uncertainty support.
  
   The rptree is stored for the table described by the nodes in
   modes[0]. All uncertainty modes stored in the modes array must
   have the same number of points and conform to the same topology
   (given by the tris array) as modes[0]. If a new table is created by
   mixing of the modes, it must be stored in modes[0] and the rptree
   regenerated using calculate_rptree.
*/
typedef struct eos_table_utri
{
  eos_table_rptree_t rptree; /* the RP tree for table lookup*/
  eos_data_points_t **modes; /* array of EOS data storage containers includes the mean values*/
  eos_table_tri_t *tris;     /* array of triangles in the table*/
  double *defaults;          /* defaults array (Z,A,rho0,T0,?)*/
  int flags[2];              /* integer flags*/
                             /* index 0: log independent variable flag*/
                             /* index 1: type of independent variables*/
  int nmodes;                /* length of the modes array */
  int np;                    /* number of points in each mode*/
  int nt;                    /* number of triangles*/
  int nbp[2];                /* number of nodes on the upper and lower boundaries*/
  int nbt[2];                /* index to first triangle for upper and lower boundaries*/
  int nbi[2];                /* index to first point for upper and lower boundaries*/
} eos_table_utri_t;

/*
   Methods to read/write from disk an eos_table_utri_t from the file "nc_name"
*/
int eostableutrireadnc(const char *nc_name,eos_table_utri_t *t,int handle);
int eostableutriwritenc(const char *nc_name,eos_table_utri_t *t);

/*
   Method to free the resources associated with a eos_table_utri_t object
*/
void eos_table_utri_delete(eos_table_utri_t *t);

/*
   Vectorized method to evaluate the table at a given set of points.
  
   input:
     t      - pointer to the EOS table data structure
     stride - size of the vectorized input being passed
     n      - number of points in the input to evaluate (starting from 0)
     x      - array of first independent variable
     y      - array of second independent variable
     type   - 1 for rho,T input vars, 2 for rho,E input vars
     scratch     - 5*stride length of scratch space
     int_scratch - stride length of scratch space
  
   output:
     f    - interpolated EOS values for all dependent variables
            structure must be allocated by the caller
     clip - length 4 array (passed by caller) containing counts of
            number of times each of the (xlow,xhigh,ylow,yhigh) clips
            were applied to keep states in the rectangular table domain
*/
int eostableutrievaluatec0(eos_table_utri_t *t,int stride,int n,double *x,double *y,eos_data_points_t *f,int type,int *clip,double *scratch,int *int_scratch);

/*
   Construct the boundary triangles for the table.
  
   The boundary triangles allow interpolation lookup via an rptree on
   a square domain when the actual domain is curved (in the second
   independent variable). State evaluation on the boundary triangles
   return values clipped to the appropriate boundary. Note, undesired
   behavior will occur if this method is call twice on the same
   table. Only the 0 mode is changed by this routine.
  
   input:
     ttype - type of table (independent vars), 2 for (rho,E)
  
   input/output:
     table - pointer to the EOS table data structure
*/
int construct_boundary_tris(eos_table_utri_t *t,int ttype);

/*
   Function in the file utri_eos_rptree.c:
  
   input:
     table - pointer to the unstructured tri eos storage structure
     ttype - flag to set independent variables for the table.
             ( 0 == rho,T ; 1 == rho,E )
  
   output:
     table->rptree - regenerated rptree based upon location of modes[0] nodes
*/
int calculate_rptree(eos_table_utri_t *table,int ttype);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
