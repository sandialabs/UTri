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
#include <math.h>
#include <string.h>

#include "eosdatapoints.h"

void eospointsfree(void *p)
{
  eos_data_points_t *points;

  points = (eos_data_points_t *)p;

  if (points != NULL) {
    points->n = 0;
    if (points->phase) free(points->phase);
    if (points->rho) free(points->rho);
    if (points->T) free(points->T);
    if (points->p) free(points->p);
    if (points->E) free(points->E);
    if (points->S) free(points->S);
    if (points->F) free(points->F);
    if (points->G) free(points->G);
    if (points->H) free(points->H);
    if (points->dpdr) free(points->dpdr);
    if (points->dpdT) free(points->dpdT);
    if (points->dEdr) free(points->dEdr);
    if (points->dEdT) free(points->dEdT);
    if (points->cs) free(points->cs);
    if (points->dcsdr) free(points->dcsdr);
    if (points->visc) free(points->visc);
    if (points->thcon) free(points->thcon);
    points->phase = NULL;
    points->rho = NULL;
    points->T = NULL;
    points->p = NULL;
    points->E = NULL;
    points->S = NULL;
    points->F = NULL;
    points->G = NULL;
    points->H = NULL;
    points->dpdr = NULL;
    points->dpdT = NULL;
    points->dEdr = NULL;
    points->dEdT = NULL;
    points->cs = NULL;
    points->dcsdr = NULL;
    points->visc = NULL;
    points->thcon = NULL;
  }

}

void *eospointsrealloc(void *p,int l)
{
  eos_data_points_t *points;
  int fail;

  points = (eos_data_points_t *)p;

  fail = 0;
  if (points == NULL) {
    if ((points = (eos_data_points_t *)malloc(sizeof(eos_data_points_t))) == NULL) fail = 1;
    if (fail) return NULL;
    
    points->n = 0;
    points->phase = NULL;
    points->rho = NULL;
    points->T = NULL;
    points->p = NULL;
    points->E = NULL;
    points->S = NULL;
    points->F = NULL;
    points->G = NULL;
    points->H = NULL;
    points->dpdr = NULL;
    points->dpdT = NULL;
    points->dEdr = NULL;
    points->dEdT = NULL;
    points->cs = NULL;
    points->dcsdr = NULL;
    points->visc = NULL;
    points->thcon = NULL;
  }

  if (points->n > 0 || l > 0) {
    points->n = l;
    if ((points->phase = (int *)realloc(points->phase,sizeof(int)*l)) == NULL) fail = 1;
    if ((points->rho = (double *)realloc(points->rho,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->T = (double *)realloc(points->T,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->p = (double *)realloc(points->p,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->E = (double *)realloc(points->E,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->S = (double *)realloc(points->S,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->F = (double *)realloc(points->F,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->G = (double *)realloc(points->G,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->H = (double *)realloc(points->H,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->dpdr = (double *)realloc(points->dpdr,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->dpdT = (double *)realloc(points->dpdT,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->dEdr = (double *)realloc(points->dEdr,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->dEdT = (double *)realloc(points->dEdT,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->cs = (double *)realloc(points->cs,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->dcsdr = (double *)realloc(points->dcsdr,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->visc = (double *)realloc(points->visc,sizeof(double)*l)) == NULL) fail = 1;
    if ((points->thcon = (double *)realloc(points->thcon,sizeof(double)*l)) == NULL) fail = 1;
  }

  if (fail) return NULL;
  else return (void *)points;
}
