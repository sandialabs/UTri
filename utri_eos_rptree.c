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

#include "utri_eos_support.h"

typedef struct node {
  double c[2]; /* location of node */
  int n;       /* number of node in table */
  int nn;      /* new location in node array */
} node_t;

typedef struct tri {
  int n[3];
  double cbl[2],cbh[2];
} tri_t;

typedef struct rtreenode {
  double split;
  int type;
  int level;
  int lt; /* index into rtreenodevec */
  int gt;
  int ntri;
  int *t;
  int *ts;
  int nn;
  int *n;
  double lb[2];
  double ub[2];
} rtreenode_t;

void gettribccoords(node_t *n,tri_t *t,double px,double py,double *l)
{
  double a[7];

  a[0] = n[t->n[0]].c[0]-n[t->n[2]].c[0];
  a[1] = n[t->n[1]].c[0]-n[t->n[2]].c[0];
  a[2] = n[t->n[2]].c[0];
  a[3] = n[t->n[0]].c[1]-n[t->n[2]].c[1];
  a[4] = n[t->n[1]].c[1]-n[t->n[2]].c[1];
  a[5] = n[t->n[2]].c[1];
  a[6] = 1./(a[0]*a[4]-a[1]*a[3]);

  l[0] = (a[1]*(a[5]-py)-(a[2]-px)*a[4])*a[6];
  l[1] = (a[3]*(a[2]-px)-(a[5]-py)*a[0])*a[6];
  l[2] = 1.-l[0]-l[1];
  if (fabs(l[0]) < 2.e-14) l[0] = 0.;
  if (fabs(l[1]) < 2.e-14) l[1] = 0.;
  if (fabs(l[2]) < 2.e-14) l[2] = 0.;
}

int rtreenode_init(rtreenode_t *n,int nt,int nn)
{
  int fail = 0;
  int i;

  if (n == NULL) {
    reportUTriError("Error in rtreenode_init: NULL input detected");
    return -1;
  }
  if (nt < 1) {
    reportUTriError("Error in rtreenode_init: tri count less than 1");
    return -1;
  }
  if (nn < 0) {
    reportUTriError("Error in rtreenode_init: node count less than 0");
    return -1;
  }

  n->ntri = nt;
  n->nn = nn;
  if ((n->t = (int *)malloc(sizeof(int)*nt)) == NULL) {
    fail = 1;
    reportUTriError("Error in rtreenode_init: failed to allocate tri list");
  }
  if ((n->ts = (int *)malloc(sizeof(int)*nt)) == NULL) {
    fail = 1;
    reportUTriError("Error in rtreenode_init: failed to allocate tri side list");
  }
  if (nn > 0) {
    if ((n->n = (int *)malloc(sizeof(int)*nn)) == NULL) {
      fail = 1;
      reportUTriError("Error in rtreenode_init: failed to allocate node list");
    }
  }
  else n->n = NULL;
  if (fail) {
    if (n->t) free(n->t);
    if (n->ts) free(n->ts);
    if (n->n) free(n->n);
    return -1;
  }

  n->lt = -1;
  n->gt = -1;
  n->type = -1;
  n->level = -1;
  n->split = 0.;
  for (i=0;i<nt;i++) {
    n->t[i] = -1;
    n->ts[i] = -1;
  }
  for (i=0;i<nn;i++) n->n[i] = -1;

  return 0;
}

void rtreenode_destroy(rtreenode_t *n)
{
  if (n == NULL) return;

  if (n->t) free(n->t);
  if (n->ts) free(n->ts);
  if (n->n) free(n->n);
}

typedef struct rtreenodevec {
  rtreenode_t *v;
  int vlen;
  int numv;
} rtreenodevec_t;

int rtreenodevec_init(rtreenodevec_t *v)
{
  if (v == NULL) {
    reportUTriError("Error in rtreenodevec_init: NULL input");
    return -1;
  }

  v->vlen = 8;
  v->numv = 0;

  if ((v->v = (rtreenode_t *)malloc(sizeof(rtreenode_t)*v->vlen)) == NULL) {
    reportUTriError("Error in rtreenodevec_init: failed to allocate memory for vector");
    return -1;
  }

  return 0;
}

void rtreenodevec_destroy(rtreenodevec_t *v)
{
  int i;
  
  if (v == NULL) return;

  for (i=0;i<v->numv;i++) {
    rtreenode_destroy(&v->v[i]);
  }
  free(v->v);
  v->v = NULL;
  v->vlen = 0;
  v->numv = 0;
   
}

int rtreenodevec_add(rtreenodevec_t *v,int *i,int nt,int nn)
{
  rtreenode_t *err;

  if (v == NULL || i == NULL) {
    reportUTriError("Error in rtreenodevec_add: NULL inputs detected");
    return -1;
  }

  /* expand vector */
  if (v->numv == v->vlen) {
    v->vlen *= 2;
    if ((err = (rtreenode_t *)realloc(v->v,sizeof(rtreenode_t)*v->vlen)) == NULL) {
      reportUTriError("Failed to reallocate rtreenode memory in rtreenodevec_add");
      return -1;
    }
    else v->v = err;
  }

  /* init node */
  if (rtreenode_init(&v->v[v->numv],nt,nn) != 0) {
    reportUTriError("Error in rtreenodevec_add: failed to initialize rtreenode");
    return -1;
  }
  *i = v->numv;
  v->numv++;

  return 0;
}
void print_tri(node_t *n,tri_t *t)
{
  printf("# n1 %d n2 %d n3 %d\n",t->n[0],t->n[1],t->n[2]);
  printf("x1 %g y1 %g\n",n[t->n[0]].c[0],n[t->n[0]].c[1]);
  printf("x2 %g y2 %g\n",n[t->n[1]].c[0],n[t->n[1]].c[1]);
  printf("x3 %g y3 %g\n",n[t->n[2]].c[0],n[t->n[2]].c[1]);
  printf("x1 %g y1 %g\n\n",n[t->n[0]].c[0],n[t->n[0]].c[1]);
}
/*

void print_tribb(tri_t *t)
{
  printf("bx1 %g by1 %g\n",t->cbl[0],t->cbl[1]);
  printf("bx2 %g by1 %g\n",t->cbh[0],t->cbl[1]);
  printf("bx3 %g by1 %g\n",t->cbh[0],t->cbh[1]);
  printf("bx4 %g by1 %g\n",t->cbl[0],t->cbh[1]);
  printf("bx1 %g by1 %g\n\n",t->cbl[0],t->cbl[1]);
}

void print_rtreenode(node_t *n,tri_t *t,rtreenode_t *rtn)
{
  int i;

  printf("\n# rtreenode\n\n");

  printf("# split bounds\n");
  printf("x1 %g y1 %g\n",rtn->lb[0],rtn->lb[1]);
  printf("x2 %g y2 %g\n",rtn->ub[0],rtn->lb[1]);
  printf("x3 %g y3 %g\n",rtn->ub[0],rtn->ub[1]);
  printf("x4 %g y4 %g\n",rtn->lb[0],rtn->ub[1]);
  printf("x1 %g y1 %g\n",rtn->lb[0],rtn->lb[1]);
  printf("\n\n");

  for (i=0;i<rtn->ntri;i++) {
    printf("# tri %d\n",i);
    print_tri(n,&t[rtn->t[i]]);
  }
  printf("\n");

  for (i=0;i<rtn->ntri;i++) {
    printf("# tri %d bb\n",i);
    print_tribb(&t[rtn->t[i]]);
  }
  printf("\n");

}
*/
int tri_side_intersections(node_t *n,tri_t *t,double *lb,double *ub,int type,double *intersections,int *ni)
{
  double x0,x1,y0,y1;
  int i;
  double s,si;
  double x,y;
  int debug = 0;

  if (debug > 0) printf("in tri_side_intersections: type %d lb %g %g ub %g %g\n",type,lb[0],lb[1],ub[0],ub[1]);

  /* side intersections */
  for (i=0;i<3;i++) {
    y1 = n[t->n[(i+1)%3]].c[1];
    y0 = n[t->n[i]].c[1];
    x1 = n[t->n[(i+1)%3]].c[0];
    x0 = n[t->n[i]].c[0];

    if (x1 == x0) { 
      /* vertical edge */
      if ( type==0 && x0 > lb[0] && x0 < ub[0]
	  && ((lb[1]-y1)*(lb[1]-y0)<0. || (ub[1]-y1)*(ub[1]-y0)<0.)) {
	if (debug > 0) printf("tsi 1 %g\n",x0);
	intersections[(*ni)++] = x0;
      }
    }
    else if (y1 == y0) {
      /* horizontal edge */
      if (type == 1 && y0 > lb[1] && y0 < ub[1]
	  && ((lb[0]-x1)*(lb[0]-x0)<0. || (ub[0]-x1)*(ub[0]-x0)<0.)) {
	if (debug > 0) printf("tsi 2 %g\n",y0);
	intersections[(*ni)++] = y0;
      }
    }
    else {
      s = (y1-y0)/(x1-x0);
      if (type == 1) {
	/* left box edge */
	y = s*(lb[0]-x0)+y0;
	if (y > lb[1]*(1.+3.e-16) && y < ub[1]*(1.-3.e-16)
            && (lb[0]-x1)*(lb[0]-x0) < 3.e-16) {
	  if (debug > 0) printf("tsi 3 %g\n",y);
	  intersections[(*ni)++] = y;
	}
	/* right box edge */
	y = s*(ub[0]-x0)+y0;
	if (y > lb[1]*(1.+3.e-16) && y < ub[1]*(1.-3.e-16)
            && (ub[0]-x1)*(ub[0]-x0) < 3.e-16) {
	  if (debug > 0) printf("tsi 4 %g\n",y);
	  intersections[(*ni)++] = y;
	}
      }
      else {
        si = 1./s; 
	/* bottom box edge */
	x = si*(lb[1]-y0)+x0;
	if (x > lb[0]*(1.+3.e-16) && x < ub[0]*(1.-3.e-16)
            && (lb[1]-y1)*(lb[1]-y0) < 3.e-16) {
	  if (debug > 0) printf("tsi 5 %g\n",x);
	  intersections[(*ni)++] = x;
	}
	/* top box edge */
	x = si*(ub[1]-y0)+x0;
	if (x > lb[0]*(1.+3.e-16) && x < ub[0]*(1.-3.e-16)
            && (ub[1]-y1)*(ub[1]-y0) < 3.e-16) {
	  if (debug > 0) printf("tsi 6 %g\n",x);
	  intersections[(*ni)++] = x;
	}
      }
    }
  }
  

  return 0;
}

int tri_rtreenode_intersection(node_t *n,tri_t *t,double *lb,double *ub,int *intersect)
{
  int i;
  int sides[3];
  int inside;
  double x,y;
  int inx,iny;
  double s;
  double x0,y0,x1,y1;
  double l[3];
  int debug = 0;

  if (debug > 0) printf("in tri_rtreenode_intersection: lb %g %g ub %g %g\n",lb[0],lb[1],ub[0],ub[1]);

  *intersect = 0;

  if (n == NULL || t == NULL || lb == NULL || ub == NULL) {
    reportUTriError("Error in tri_rtreenode_intersection: NULL inputs detected");
    return -1;
  }

  /* check location of triangle vertices in relation to square */
  inside = 0;
  sides[0] = -1;
  sides[1] = -1;
  sides[2] = -1;
  for (i=0;i<3;i++) {
    x = n[t->n[i]].c[0];
    y = n[t->n[i]].c[1];
    if (x > lb[0] && x < ub[0]) inx = 1;
    else inx = 0;
    if (y > lb[1] && y < ub[1]) iny = 1;
    else iny = 0;

    if (inx && iny) inside++;
    else if (y == lb[1]) {
      if (inx) sides[i] = 0;
      else if (x == lb[0]) sides[i] = 4;
      else if (x == ub[0]) sides[i] = 5;
    }
    else if (y == ub[1]) {
      if (inx) sides[i] = 2;
      else if (x == ub[0]) sides[i] = 6;
      else if (x == lb[0]) sides[i] = 7;
    }
    else if (x == ub[0]) {
      if (iny) sides[i] = 1;
    }
    else if (x == lb[0]) {
      if (iny) sides[i] = 3;
    }
  }
  /* triangle point inside */
  if (inside > 0) {
    *intersect = 1;
    return 0;
  }
  if (debug > 0) for (i=0;i<3;i++) printf("side %d %d\n",i,sides[i]);
  /* triangle inscribed in square */
  if (sides[0] != -1 && sides[1] != -1 && sides[2] != -1) {
    *intersect = 2;
    return 0;
  }
  for (i=0;i<3;i++) {
    /* triangle point along a side */
    if (sides[i] == 0 && 
	(n[t->n[(i+1)%3]].c[1] > lb[1] || n[t->n[(i+2)%3]].c[1] > lb[1])) {
      *intersect = 3;
      return 0;
    }
    if (sides[i] == 1 && 
	(n[t->n[(i+1)%3]].c[0] < ub[0] || n[t->n[(i+2)%3]].c[0] < ub[0])) {
      *intersect = 4;
      return 0;
    }
    if (sides[i] == 2 && 
	(n[t->n[(i+1)%3]].c[1] < ub[1] || n[t->n[(i+2)%3]].c[1] < ub[1])) {
      *intersect = 5;
      return 0;
    }
    if (sides[i] == 3 && 
	(n[t->n[(i+1)%3]].c[0] > lb[0] || n[t->n[(i+2)%3]].c[0] > lb[0])) {
      *intersect = 6;
      return 0;
    }
    /* triangle point in a corner */
    if (sides[i] == 4 &&
	((n[t->n[(i+2)%3]].c[1] > lb[1] && n[t->n[(i+1)%3]].c[0] > lb[0]) ||
	 (n[t->n[(i+1)%3]].c[1] > lb[1] && n[t->n[(i+1)%3]].c[0] > lb[0]) ||
	 (n[t->n[(i+2)%3]].c[1] > lb[1] && n[t->n[(i+2)%3]].c[0] > lb[0]))) {
      *intersect = 7;
      return 0;
    }
    /* triangle point in a corner */
    if (sides[i] == 5 &&
	((n[t->n[(i+1)%3]].c[1] > lb[1] && n[t->n[(i+2)%3]].c[0] < ub[0]) ||
	 (n[t->n[(i+1)%3]].c[1] > lb[1] && n[t->n[(i+1)%3]].c[0] < ub[0]) ||
	 (n[t->n[(i+2)%3]].c[1] > lb[1] && n[t->n[(i+2)%3]].c[0] < ub[0]))) {
      *intersect = 8;
      return 0;
    }
    /* triangle point in a corner */
    if (sides[i] == 6 &&
	((n[t->n[(i+2)%3]].c[1] < ub[1] && n[t->n[(i+1)%3]].c[0] < ub[0]) ||
	 (n[t->n[(i+1)%3]].c[1] < ub[1] && n[t->n[(i+1)%3]].c[0] < ub[0]) ||
	 (n[t->n[(i+2)%3]].c[1] < ub[1] && n[t->n[(i+2)%3]].c[0] < ub[0]))) {
      *intersect = 9;
      return 0;
    }
    /* triangle point in a corner */
    if (sides[i] == 7 &&
	((n[t->n[(i+1)%3]].c[1] < ub[1] && n[t->n[(i+2)%3]].c[0] > lb[0]) ||
	 (n[t->n[(i+1)%3]].c[1] < ub[1] && n[t->n[(i+1)%3]].c[0] > lb[0]) ||
	 (n[t->n[(i+2)%3]].c[1] < ub[1] && n[t->n[(i+2)%3]].c[0] > lb[0]))) {
      *intersect = 10;
      return 0;
    }
  }
  /* triangle with point on side failing above tests is outside */
  if (sides[0] != -1 || sides[1] != -1 || sides[2] != -1) return 0;
  if (debug > 0) printf("checking side intersections\n");
  /* side intersections */
  for (i=0;i<3;i++) {
    y1 = n[t->n[(i+1)%3]].c[1];
    y0 = n[t->n[i]].c[1];
    x1 = n[t->n[(i+1)%3]].c[0];
    x0 = n[t->n[i]].c[0];
    if (debug > 0) printf("triside p1 %g %g p2 %g %g\n",x0,y0,x1,y1);

    if (x1 == x0) {
      /* vertical edge */
      if (x0 > lb[0] && x0 < ub[0] 
	  && ((lb[1]-y1)*(lb[1]-y0)<0. || (ub[1]-y1)*(ub[1]-y0)<0.)) {
	*intersect = 11;
	return 0;
      }
    }
    else if (y1 == y0) {
      /* horizontal edge */
      if (y0 > lb[1] && y0 < ub[1]
	  && ((lb[0]-x1)*(lb[0]-x0)<0. || (ub[0]-x1)*(ub[0]-x0)<0.)) {
	*intersect = 12;
	return 0;
      }
    }
    else {
      s = (y1-y0)/(x1-x0);
      /* left box edge */
      y = s*(lb[0]-x0)+y0;
      if (debug > 0) printf("lbe y %g\n",y);
      if (y >= lb[1] && y <= ub[1]
	  && (lb[0]-x1)*(lb[0]-x0) < 0.) {
	*intersect = 13;
	return 0;
      }
      /* right box edge */
      y = s*(ub[0]-x0)+y0;
      if (debug > 0) printf("rbe y %g\n",y);
      if (y >= lb[1] && y <= ub[1]
	  && (ub[0]-x1)*(ub[0]-x0) < 0.) {
	*intersect = 14;
	return 0;
      }
      /* bottom box edge */
      x = (1./s)*(lb[1]-y0)+x0;
      if (debug > 0) printf("bbe x %g\n",x);
      if (x >= lb[0] && x <= ub[0]
	  && (lb[1]-y1)*(lb[1]-y0) < 0.) {
	*intersect = 15;
	return 0;
      }
/*
         top box edge
         this test is redundant as a triangle will intersect any two sides
         all 6 cases of two side intersection are covered by the above three tests
        x = (1./s)*(ub[1]-y0)+x0;
        if (x >= lb[0] && x <= ub[0]
            && (ub[1]-y1)*(ub[1]-y0) < 0.) {
          *intersect = 16;
          return 0;
        }
*/
    }
  }
  /* square inside triangle */
  inside = 0;
  gettribccoords(n,t,lb[0],lb[1],l);
  if ((l[0] >= 0. && l[0] <= 1.) &&
      (l[1] >= 0. && l[1] <= 1.) &&
      (l[2] >= 0. && l[2] <= 1.)) inside++;    
  gettribccoords(n,t,ub[0],lb[1],l);
  if ((l[0] >= 0. && l[0] <= 1.) &&
      (l[1] >= 0. && l[1] <= 1.) &&
      (l[2] >= 0. && l[2] <= 1.)) inside++;    
  gettribccoords(n,t,ub[0],ub[1],l);
  if ((l[0] >= 0. && l[0] <= 1.) &&
      (l[1] >= 0. && l[1] <= 1.) &&
      (l[2] >= 0. && l[2] <= 1.)) inside++;    
  gettribccoords(n,t,lb[0],ub[1],l);
  if ((l[0] >= 0. && l[0] <= 1.) &&
      (l[1] >= 0. && l[1] <= 1.) &&
      (l[2] >= 0. && l[2] <= 1.)) inside++;    
  if (inside == 4) {
    *intersect = 17;
    return 0;
  }

  return 0;
}

node_t *ds_n;
int ds_type;
int ds_tsort(int *a,int *b)
{
  if (ds_n[*a].c[ds_type] > ds_n[*b].c[ds_type]) return 1;
  else if (ds_n[*a].c[ds_type] < ds_n[*b].c[ds_type]) return -1;
  else return 0;
}

int dsort(double *a,double *b)
{
  if (*a > *b) return 1;
  else if (*a < *b) return -1;
  else return 0;
}

int dosplit(node_t *n,tri_t *t,rtreenodevec_t *r,int ri)
{
  int i,j;
  int totl,toth;
  rtreenode_t *rtn;
  rtreenode_t *rtnlt;
  rtreenode_t *rtngt;
  int type;
  int lt,gt;
  int intersect;
  double lb1[2],ub1[2];
  double lb2[2],ub2[2];
  double splits[3];
  int totls[3],toths[3],startns[3];
  int startn;
  int totnl,totnh;
  int insidenodes;
  double *splitlocs;
  int nsl;
  int tot;
  /* left and right inside node */
  int inl,inr;
  int inm;
  int debug = 0;

  /* node to split */
  rtn = &r->v[ri];

  /* type of split */
  type = rtn->type;

  if (debug > 0) {
    printf("# dosplit starting on %d nodes\n",rtn->nn);
    printf("# bounding box:\n%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n\n",rtn->lb[0],rtn->lb[1],
	   rtn->ub[0],rtn->lb[1],rtn->ub[0],rtn->ub[1],rtn->lb[0],rtn->ub[1],rtn->lb[0],rtn->lb[1]);
  }

  /* setup sort routines */
  ds_n = n;
  ds_type = type;

  /* setup bounds */
  lb1[0] = rtn->lb[0];
  lb1[1] = rtn->lb[1];
  ub1[0] = rtn->ub[0];
  ub1[1] = rtn->ub[1];
  lb2[0] = rtn->lb[0];
  lb2[1] = rtn->lb[1];
  ub2[0] = rtn->ub[0];
  ub2[1] = rtn->ub[1];

  /* sort nodes according to type */
  qsort(rtn->n,rtn->nn,sizeof(int),(int (*)(const void *,const void *))ds_tsort);

  /* check if nodes allow standard split technique */
  insidenodes = 0;
  inl = -1;
  inr = -1;
  for (i=0;i<rtn->nn;i++) {
    if (n[rtn->n[i]].c[type] > rtn->lb[type] && n[rtn->n[i]].c[type] < rtn->ub[type]) {
      insidenodes = 1;
      if (inl < 0) inl = i;
      inr = i;
    }
    else if (insidenodes == 1) break;
  }
  if (debug > 0) printf("# insidenodes %d\n",insidenodes);
  
  if (insidenodes == 0) {
    /* upper bound to number of locations */
    nsl = rtn->ntri*6; 
    if ((splitlocs = (double *)malloc(sizeof(double)*nsl)) == NULL) {
      reportUTriError("Failed to allocate memory for split locations in dosplit");
      return -1;
    }
    nsl = 0;

    /* find intersections of triangles with splitting sides */
    for (i=0;i<rtn->ntri;i++) {
      if (tri_side_intersections(n,&t[rtn->t[i]],rtn->lb,rtn->ub,type,splitlocs,&nsl) != 0) {
	reportUTriError("tri_side_intersections failed");
	return -1;
      }
    }
    
    if (debug > 0) {
      for (i=0;i<nsl;i++) {
	printf("# sl %d type %d val %g\n",i,type,splitlocs[i]);
      }
    }
    
    /* sort intersections */
    qsort(splitlocs,nsl,sizeof(double),(int (*)(const void *,const void *))dsort);

    if (debug > 0) {
      printf("# sorted sl\n");
      for (i=0;i<nsl;i++) {
	printf("# sl %d type %d val %g\n",i,type,splitlocs[i]);
      }
    }

    /* prune off any non-pair intersections */
    j=0;
    tot=1;
    for (i=1;i<nsl;i++) {
      if (fabs(splitlocs[i]-splitlocs[j]) < 1.e-14*fabs(splitlocs[i])) {
	if (tot == 1) {
	  splitlocs[++j] = splitlocs[i];
	  tot = 0;
	}
      }
      else {
	if (tot == 1) {
	  splitlocs[j] = splitlocs[i];
	}
	else {
	  splitlocs[++j] = splitlocs[i];
	  tot = 1;
	}
      }
    }
    nsl = j+1-tot;
    if (nsl < 2) nsl = 0;

    /* prune off intersections on edge of boundary */
    j=0;
    for (i=0;i<nsl;i++) {
      if (debug > 0) printf("#sl %g lb %g ub %g lbd %g ubd %g\n",splitlocs[i],rtn->lb[type],rtn->ub[type],splitlocs[i]-rtn->lb[type],splitlocs[i]-rtn->ub[type]);
      if ((splitlocs[i]-rtn->lb[type])>1.e-14*(splitlocs[i]+rtn->lb[type])
          && (splitlocs[i]-rtn->ub[type])>1.e-14*(splitlocs[i]+rtn->ub[type])) {
	splitlocs[j++] = splitlocs[i];
      }
    }
    nsl = j;

    if (debug > 0) {
      printf("# culled sl lb %g ub %g\n",rtn->lb[type],rtn->ub[type]);
      for (i=0;i<nsl;i++) {
	printf("# sl %d type %d val %g\n",i,type,splitlocs[i]);
      }
    }

    splits[0] = 0.;
    totls[0] = -1;
    toths[0] = -1;
    /* compute best split */
    if (nsl == 0) {
      splits[0] = (rtn->lb[type]+rtn->ub[type])/2.;
      totls[0] = rtn->ntri;
      toths[0] = rtn->ntri;
    }
    else {
      for (j=0;j<nsl/2+1;j++) {

	if (j==0) rtn->split = (rtn->lb[type]+4.*splitlocs[0])/5.;
	else if (j==nsl/2) rtn->split = (rtn->ub[type]+4.*splitlocs[nsl-1])/5.;
	else rtn->split = (splitlocs[2*j-1]+splitlocs[2*j])/2.;

	if (debug > 0) printf("# j %d nsl %d split %g\n",j,nsl,rtn->split);
	ub1[type] = rtn->split;
	lb2[type] = rtn->split;
      
	/* decide which tris are on which side */
	totl = 0;
	toth = 0;
	for (i=0;i<rtn->ntri;i++) {
	  if (debug > 0) print_tri(n,&t[rtn->t[i]]);
	  if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb1,ub1,&intersect) != 0) {
	    reportUTriError("tri_rtreenode_intersection failed");
	    return -1;
	  }
	  if (debug > 0) printf("# intersect %d\n",intersect);
	  if (intersect > 0) {
	    rtn->ts[i] = 0;
	    totl++;
	  }
	  else rtn->ts[i] = -1;
	  
	  if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb2,ub2,&intersect) != 0) {
	    reportUTriError("tri_rtreenode_intersection failed");
	    return -1;
	  }
	  if (debug > 0) printf("# intersect %d\n",intersect);
	  if (intersect > 0) {
	    if (rtn->ts[i] == 0) rtn->ts[i] = 2;
	    else rtn->ts[i] = 1;
	    toth++;
	  }
	  if (debug > 0) printf("# i %d ts %d\n",i,rtn->ts[i]);
	}
	if (debug > 0) printf("# totl %d toth %d split %g\n",totl,toth,rtn->split);
	
	if (j==0 ||
	    (totl+toth < totls[0]+toths[0] && fabs(toth-totl) <= fabs(totls[0]-toths[0])) ||
	    (totl+toth <= totls[0]+toths[0] && fabs(toth-totl) < fabs(totls[0]-toths[0]))) {
	  /* record best split info */
	  splits[0] = rtn->split;
	  totls[0] = totl;
	  toths[0] = toth;
	}
      }
    }

    /* if no improvement, default to halving interval */
    if (totls[0] == toths[0] && totls[0] == rtn->ntri) {
      splits[0] = (rtn->lb[type]+rtn->ub[type])/2.;
    }

    /* recalculate ts array for chosen split */
    rtn->split = splits[0];

    ub1[type] = rtn->split;
    lb2[type] = rtn->split;

    /* decide which tris are on which side */
    totl = 0;
    toth = 0;
    for (i=0;i<rtn->ntri;i++) {
      if (debug > 0) print_tri(n,&t[rtn->t[i]]);
      if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb1,ub1,&intersect) != 0) {
	reportUTriError("tri_rtreenode_intersection failed");
	return -1;
      }
      if (intersect > 0) {
	rtn->ts[i] = 0;
	totl++;
      }
      else rtn->ts[i] = -1;
      
      if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb2,ub2,&intersect) != 0) {
	reportUTriError("tri_rtreenode_intersection failed");
	return -1;
      }
      if (intersect > 0) {
	if (rtn->ts[i] == 0) rtn->ts[i] = 2;
	else rtn->ts[i] = 1;
	toth++;
      }
    }
    if (debug > 0) printf("# totl %d toth %d split %g\n",totl,toth,rtn->split);

    /* count number of nodes */
    totnl = 0;
    totnh = 0;
    for (i=0;i<rtn->nn;i++) {
      if (n[rtn->n[i]].c[type] < rtn->split) totnl++;
      else totnh++;
    }

    if (debug > 0) {
      printf("# totnl %d totnh %d\n",totnl,totnh);

      if (totl == 0 || toth == 0) {
	printf("# type %d bounds %g %g split %g\n",type,rtn->lb[type],rtn->ub[type],rtn->split);
	printf("# rtn bounds x %g %g y %g %g\n",rtn->lb[0],rtn->ub[0],rtn->lb[1],rtn->ub[1]);
      }
    }

    free(splitlocs);
  }
  else {

    if (debug > 0) {
      printf("# inl %d inr %d\n",inl,inr);
      for (i=0;i<rtn->nn;i++) printf("# in %d xy %g\n",i,n[rtn->n[i]].c[type]);
    }

    /* split at midpoint of inside nodes */
    inm = inl+inr+1;
    if (inm % 2) {
      rtn->split = n[rtn->n[inm/2]].c[type];
    }
    else {
      rtn->split = (n[rtn->n[inm/2-1]].c[type]+
		    n[rtn->n[inm/2]].c[type])/2.;
    }

    ub1[type] = rtn->split;
    lb2[type] = rtn->split;

    /* decide which tris are on which side */
    totl = 0;
    toth = 0;
    for (i=0;i<rtn->ntri;i++) {
      if (debug > 0) print_tri(n,&t[rtn->t[i]]);
      if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb1,ub1,&intersect) != 0) {
	reportUTriError("tri_rtreenode_intersection failed");
	return -1;
      }
      if (intersect > 0) {
	rtn->ts[i] = 0;
	totl++;
      }
      else rtn->ts[i] = -1;
      
      if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb2,ub2,&intersect) != 0) {
	reportUTriError("tri_rtreenode_intersection failed");
	return -1;
      }
      if (intersect > 0) {
	if (rtn->ts[i] == 0) rtn->ts[i] = 2;
	else rtn->ts[i] = 1;
	toth++;
      }
    }
    if (debug > 0) printf("# totl %d toth %d\n",totl,toth);
    
    /* record initial split info */
    for (i=0;i<3;i++) {
      splits[i] = rtn->split;
      totls[i] = totl;
      toths[i] = toth;
      if (inm%2) startns[i] = inm/2;
      else startns[i] = -1;
    }

    if (debug > 0) printf("# initsplit %g totl %d toth %d startns %d\n",splits[0],totls[0],toths[0],startns[0]);
    
    /* try reducing the split, choosing node locations */
    startn = inm/2;
    while (1) {
      startn--;
      if (startn < inl) break;
      
      rtn->split = n[rtn->n[startn]].c[type];
      if (debug > 0) printf("# split %g\n",rtn->split);
      ub1[type] = rtn->split;
      lb2[type] = rtn->split;
      
      /* decide which tris are on which side */
      totl = 0;
      toth = 0;
      for (i=0;i<rtn->ntri;i++) {
	if (debug > 0) print_tri(n,&t[rtn->t[i]]);
	if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb1,ub1,&intersect) != 0) {
	  reportUTriError("tri_rtreenode_intersection failed");
	  return -1;
	}
	if (intersect > 0) {
	  rtn->ts[i] = 0;
	  totl++;
	}
	else rtn->ts[i] = -1;
	
	if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb2,ub2,&intersect) != 0) {
	  reportUTriError("tri_rtreenode_intersection failed");
	  return -1;
	}
	if (intersect > 0) {
	  if (rtn->ts[i] == 0) rtn->ts[i] = 2;
	  else rtn->ts[i] = 1;
	  toth++;
	}
	if (debug > 0) printf("# i %d ts %d\n",i,rtn->ts[i]);
      }
      if (debug > 0) printf("# totl %d toth %d startn %d\n",totl,toth,startn);
      
      if ((totl+toth < totls[1]+toths[1] && fabs(toth-totl) <= fabs(totls[1]-toths[1])) ||
	  (totl+toth <= totls[1]+toths[1] && fabs(toth-totl) < fabs(totls[1]-toths[1]))) {
	/* record best split info */
	splits[1] = rtn->split;
	totls[1] = totl;
	toths[1] = toth;
	startns[1] = startn;
      }
      else {
	break;
      }
    }
    
    /* try increasing the split, choosing node locations */
    startn = inm/2-1+(inm%2);
    while (1) {
      startn++;
      if (startn > inr) break;
      
      rtn->split = n[rtn->n[startn]].c[type];
      
      ub1[type] = rtn->split;
      lb2[type] = rtn->split;
      
      /* decide which tris are on which side */
      totl = 0;
      toth = 0;
      for (i=0;i<rtn->ntri;i++) {
	if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb1,ub1,&intersect) != 0) {
	  reportUTriError("tri_rtreenode_intersection failed");
	  return -1;
	}
	if (intersect > 0) {
	  rtn->ts[i] = 0;
	  totl++;
	}
	else rtn->ts[i] = -1;
	
	if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb2,ub2,&intersect) != 0) {
	  reportUTriError("tri_rtreenode_intersection failed");
	  return -1;
	}
	if (intersect > 0) {
	  if (rtn->ts[i] == 0) rtn->ts[i] = 2;
	  else rtn->ts[i] = 1;
	  toth++;
	}
	if (debug > 0) printf("# i %d ts %d intersect %d\n",i,rtn->ts[i],intersect);
      }
      if (debug > 0) printf("# totl %d toth %d startn %d\n",totl,toth,startn);
      
      if ((totl+toth < totls[2]+toths[2] && fabs(toth-totl) <= fabs(totls[2]-toths[2])) ||
	  (totl+toth <= totls[2]+toths[2] && fabs(toth-totl) < fabs(totls[2]-toths[2]))) {
	/* record best split info */
	splits[2] = rtn->split;
	totls[2] = totl;
	toths[2] = toth;
	startns[2] = startn;
      }
      else {
	break;
      }
    }
    
    /* choose best split */
    if ((totls[2]+toths[2] < totls[0]+toths[0] && fabs(toths[2]-totls[2]) <= fabs(totls[0]-toths[0])) ||
	(totls[2]+toths[2] <= totls[0]+toths[0] && fabs(toths[2]-totls[2]) < fabs(totls[0]-toths[0]))) {
      /* record best split info */
      splits[0] = splits[2];
      totls[0] = totls[2];
      toths[0] = toths[2];
      startns[0] = startns[2];
    }
    
    /* choose best split */
    if ((totls[1]+toths[1] < totls[0]+toths[0] && fabs(toths[1]-totls[1]) <= fabs(totls[0]-toths[0])) ||
	(totls[1]+toths[1] <= totls[0]+toths[0] && fabs(toths[1]-totls[1]) < fabs(totls[0]-toths[0]))) {
      /* record best split info */
      splits[0] = splits[1];
      totls[0] = totls[1];
      toths[0] = toths[1];
      startns[0] = startns[1];
    }
    
    /* recalculate ts array for chosen split */
    rtn->split = splits[0];

    ub1[type] = rtn->split;
    lb2[type] = rtn->split;

    /* decide which tris are on which side */
    totl = 0;
    toth = 0;
    for (i=0;i<rtn->ntri;i++) {
      if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb1,ub1,&intersect) != 0) {
	reportUTriError("tri_rtreenode_intersection failed");
	return -1;
      }
      if (intersect > 0) {
	rtn->ts[i] = 0;
	totl++;
      }
      else rtn->ts[i] = -1;
      if (debug > 0) printf("# i %d ts %d intersect %d\n",i,rtn->ts[i],intersect);
      
      if (tri_rtreenode_intersection(n,&t[rtn->t[i]],lb2,ub2,&intersect) != 0) {
	reportUTriError("tri_rtreenode_intersection failed");
	return -1;
      }
      if (intersect > 0) {
	if (rtn->ts[i] == 0) rtn->ts[i] = 2;
	else rtn->ts[i] = 1;
	toth++;
      }
      if (debug > 0) printf("# i %d ts %d intersect %d\n",i,rtn->ts[i],intersect);
    }
    if (debug > 0) printf("# totl %d toth %d splitn %d split %g\n",totl,toth,startns[0],splits[0]);

    if (startns[0] == -1) {
      /* split between inm/2 and inm/2-1 nodes */
      totnl = inm/2;
      totnh = rtn->nn-inm/2;
    }
    else {
      totnl = startns[0]+1;
      totnh = rtn->nn-startns[0];
    }
  }

  /* if multiple nodes at split location, add them to the split */
  if (debug > 0) printf("# a nn %d totnl %d totnh %d\n",rtn->nn,totnl,totnh);
  while (totnl < rtn->nn) {
    if (n[rtn->n[totnl]].c[type] == rtn->split) totnl++;
    else break;
  }
  while (totnh < rtn->nn) {
    if (n[rtn->n[rtn->nn-1-totnh]].c[type] == rtn->split) totnh++;
    else break;
  }
  if (debug > 0) printf("# b nn %d totnl %d totnh %d\n",rtn->nn,totnl,totnh);

  /* make split */
  if (rtreenodevec_add(r,&lt,totl,totnl) != 0) {
    reportUTriError("Error in dosplit: failed to add lt node to vector");
    return -1;
  }
  if (rtreenodevec_add(r,&gt,toth,totnh) != 0) {
    reportUTriError("Error in dosplit: failed to add gt node to vector");
    return -1;
  }

  /* reinitialize rtn in case nodes became reallocated */
  rtn = &r->v[ri];

  /* save split locations */
  rtn->lt = lt;
  rtn->gt = gt;
  rtnlt = &r->v[rtn->lt];
  rtngt = &r->v[rtn->gt];

  /* update bounds */
  rtnlt->lb[0] = rtn->lb[0];
  rtnlt->lb[1] = rtn->lb[1];
  rtnlt->ub[0] = rtn->ub[0];
  rtnlt->ub[1] = rtn->ub[1];
  rtngt->lb[0] = rtn->lb[0];
  rtngt->lb[1] = rtn->lb[1];
  rtngt->ub[0] = rtn->ub[0];
  rtngt->ub[1] = rtn->ub[1];
  rtnlt->ub[type] = rtn->split;
  rtngt->lb[type] = rtn->split;

  /* populate splits */
  totl = 0;
  toth = 0;
  if (debug > 0) printf("# ns %d %d\n",rtnlt->ntri,rtngt->ntri);
  for (i=0;i<rtn->ntri;i++) {
    if (rtn->ts[i] == 0) rtnlt->t[totl++] = rtn->t[i];
    else if (rtn->ts[i] == 1) rtngt->t[toth++] = rtn->t[i];
    else if (rtn->ts[i] == 2) {
      rtnlt->t[totl++] = rtn->t[i];
      rtngt->t[toth++] = rtn->t[i];
    }
  }
  rtnlt->ntri = totl;
  rtngt->ntri = toth;
  if (debug > 0) {
    printf("# na %d %d\n",rtnlt->ntri,rtngt->ntri);
    printf("# nn %d totnl %d totnh %d\n",rtn->nn,totnl,totnh);
  }
  for (i=0;i<totnl;i++) rtnlt->n[i] = rtn->n[i];
  for (i=0;i<totnh;i++) rtngt->n[i] = rtn->n[rtn->nn-1-i];
  rtnlt->type = (type+1)%2;
  rtngt->type = (type+1)%2;
  rtnlt->level = rtn->level+1;
  rtngt->level = rtn->level+1;

  if (debug > 0) printf("\n\n");

  return 0;
}

void setbb(tri_t *t,node_t *n)
{
  t->cbl[0] = n[t->n[0]].c[0];
  t->cbl[1] = n[t->n[0]].c[1];
  t->cbh[0] = n[t->n[0]].c[0];
  t->cbh[1] = n[t->n[0]].c[1];
  if (n[t->n[1]].c[0] < t->cbl[0]) t->cbl[0] = n[t->n[1]].c[0];
  if (n[t->n[1]].c[0] > t->cbh[0]) t->cbh[0] = n[t->n[1]].c[0];
  if (n[t->n[1]].c[1] < t->cbl[1]) t->cbl[1] = n[t->n[1]].c[1];
  if (n[t->n[1]].c[1] > t->cbh[1]) t->cbh[1] = n[t->n[1]].c[1];
  if (n[t->n[2]].c[0] < t->cbl[0]) t->cbl[0] = n[t->n[2]].c[0];
  if (n[t->n[2]].c[0] > t->cbh[0]) t->cbh[0] = n[t->n[2]].c[0];
  if (n[t->n[2]].c[1] < t->cbl[1]) t->cbl[1] = n[t->n[2]].c[1];
  if (n[t->n[2]].c[1] > t->cbh[1]) t->cbh[1] = n[t->n[2]].c[1];
}

typedef struct intqnode
{
  int i;
  struct intqnode *next;
} intqnode_t;

typedef struct intq
{
  intqnode_t *front;
  intqnode_t *back;
  int tot;
} intq_t;

int intq_init(intq_t *q)
{
  if (q == NULL) {
    reportUTriError("Error in intq_init: NULL input detected");
    return -1;
  }

  q->front = NULL;
  q->back = NULL;
  q->tot = 0;

  return 0;
}

void intq_destroy(intq_t *q)
{
  intqnode_t *cur;

  if (q == NULL) return;

  while (q->front != NULL) {
    cur = q->front;
    q->front = cur->next;
    if (q->back == cur) q->back = q->front;
    free(cur);
  }    
  q->tot = 0;
}

int intq_pushback(intq_t *q,int i)
{
  intqnode_t *new;

  if (q == NULL) {
    reportUTriError("Error in intq_pushback: NULL input detected");
    return -1;
  }
  
  if ((new = (intqnode_t *)malloc(sizeof(intqnode_t))) == NULL) {
    reportUTriError("Error in intq_pushback: failed to allocate node memory");
    return -1;
  }

  new->i = i;
  new->next = NULL;
  if (q->back == NULL) {
    q->front = new;
    q->back = new;
    q->tot = 1;
  }
  else {
    q->back->next = new;
    q->back = new;
    q->tot++;
  }

  return 0;
}

int intq_popfront(intq_t *q,int *i)
{
  intqnode_t *old;

  if (q == NULL) {
    reportUTriError("Error in intq_popfront: NULL input detected");
    return -1;
  }
  if (q->front == NULL) {
    reportUTriError("Error in intq_popfront: queue empty");
    return -1;
  }
  
  old = q->front;
  q->front = old->next;
  if (old == q->back) q->back = q->front;
  q->tot--;

  *i = old->i;
  free(old);

  return 0;
  
}

int save_utri_rtree(eos_table_utri_t *t,tri_t *tris,int ntris,node_t *nodes,int nnodes,rtreenodevec_t *rtnv)
{
  int i,j,k;
  tri_t *x;
  eos_table_tri_t *y;
  eos_table_rptree_t *r;
  int tot;
  int tricount;
  double invdet;
  char mesg[80];

  /* update the precomputed barycentric coordinate constants */
  for (i=0;i<ntris;i++) {
    x = &tris[i];
    y = &t->tris[i];
    y->a[0] = nodes[x->n[0]].c[0]-nodes[x->n[2]].c[0];
    y->a[1] = nodes[x->n[1]].c[0]-nodes[x->n[2]].c[0];
    y->a[2] = nodes[x->n[2]].c[0];
    y->a[3] = nodes[x->n[0]].c[1]-nodes[x->n[2]].c[1];
    y->a[4] = nodes[x->n[1]].c[1]-nodes[x->n[2]].c[1];
    y->a[5] = nodes[x->n[2]].c[1];
    invdet = 1./(y->a[0]*y->a[4]-y->a[1]*y->a[3]);
    y->a[0] *= invdet;
    y->a[1] *= invdet;
    y->a[3] *= invdet;
    y->a[4] *= invdet;
  }

  /* build rptree */
  r = &t->rptree;

  /* bounds */
  r->xbounds[0] = rtnv->v[0].lb[0];
  r->xbounds[1] = rtnv->v[0].ub[0];
  r->ybounds[0] = rtnv->v[0].lb[1];
  r->ybounds[1] = rtnv->v[0].ub[1];

  /* levels number */
  tot = 2;
  r->slevels = 0;
  while (tot < rtnv->numv+1) {
    tot *= 2;
    r->slevels++;
  }
  if (tot != rtnv->numv+1) {
    reportUTriError("bad number of rtree nodes in save_utri_rtree");
    return -1;
  }

  /* splits */
  if (r->splits) free(r->splits);
  if ((r->splits = (double *)malloc(sizeof(double)*(tot/2-1))) == NULL) {
    reportUTriError("failed to allocate splits memory in save_utri_rtree");
    return -1;
  }
  for (i=0;i<tot/2-1;i++) r->splits[i] = rtnv->v[i].split;

  /* tri index */
  if (r->trilookup) free(r->trilookup);
  if ((r->trilookup = (int *)malloc(sizeof(int)*tot)) == NULL) {
    reportUTriError("failed to allocate trilookup memory in save_utri_rtree");
    return -1;
  }
  tricount = 0;
  for (i=0;i<tot/2;i++) {
    r->trilookup[2*i+0] = rtnv->v[tot/2-1+i].ntri;
    r->trilookup[2*i+1] = tricount;
    tricount += r->trilookup[2*i+0];
  }

  /* tri list */
  if (r->trilist) free(r->trilist);
  if ((r->trilist = (int *)malloc(sizeof(int)*tricount)) == NULL) {
    reportUTriError("failed to allocate trilist memory in save_utri_rtree");
    return -1;
  }
  k = 0;
  for (i=0;i<tot/2;i++) {
    for (j=0;j<rtnv->v[tot/2-1+i].ntri;j++) {
      r->trilist[k++] = rtnv->v[tot/2-1+i].t[j];
    }
  }
  if (k != tricount) {
    snprintf(mesg,80,"error filling trilist %d %d in save_utri_rtree",k,tricount);
    reportUTriError(mesg);
    return -1;
  }

  return 0;
}

int create_utri_rtree(rtreenodevec_t *rtnv,node_t *n,int nn,tri_t *tris,int nt)
{
  int i;
  int rti;
  rtreenode_t *rt;
  rtreenode_t *rtlt;
  /*rtreenode_t *rtgt; */
  intq_t nodeq;
  int nextnode;
  int maxlevel;
  /*printf("nn %d nt %d\n",nn,nt); */

  maxlevel = ceil(log(nn)/log(2.))+2;

  if (rtreenodevec_init(rtnv) != 0) {
    reportUTriError("Error: rtreenodevec_init failed in create_utri_rtree");
    return -1;
  }

  if (rtreenodevec_add(rtnv,&rti,nt,nn) != 0) {
    reportUTriError("Error: rtreenodevec_add failed in create_utri_rtree");
    return -1;
  }

  rt = &rtnv->v[rti];
  for (i=0;i<nn;i++) rt->n[i] = i;
  for (i=0;i<nt;i++) rt->t[i] = i;
  rt->type = 0;
  rt->level = 0;
  /*
     determine initial bounding box from node bounds
  */
  rt->lb[0] = 9.99e99;
  rt->lb[1] = 9.99e99;
  rt->ub[0] = -9.99e99;
  rt->ub[1] = -9.99e99;
  for (i=0;i<nn;i++) {
    if (n[i].c[0] < rt->lb[0]) rt->lb[0] = n[i].c[0];
    if (n[i].c[1] < rt->lb[1]) rt->lb[1] = n[i].c[1];
    if (n[i].c[0] > rt->ub[0]) rt->ub[0] = n[i].c[0];
    if (n[i].c[1] > rt->ub[1]) rt->ub[1] = n[i].c[1];
  }

  if (intq_init(&nodeq) != 0) {
    reportUTriError("Error: intq_init failed in create_utri_rtree");
    return -1;
  }
  if (intq_pushback(&nodeq,rti) != 0) {
    reportUTriError("Error: intq_pushback failed in create_utri_rtree");
    return -1;
  }
  
  while (nodeq.tot > 0) {
    if (intq_popfront(&nodeq,&nextnode) != 0) {
      reportUTriError("Error: intq_popfront failed in create_utri_rtree");
      return -1;
    }

    if (dosplit(n,tris,rtnv,nextnode) != 0) {
      reportUTriError("Error: dosplit failed in create_utri_rtree");
      return -1;
    }
    
    rt = &rtnv->v[nextnode];
    rtlt = &rtnv->v[rt->lt];
    /*
    rtgt = &rtnv->v[rt->gt];
    printf("# split %g level %d type %d nn %d %d nt %d %d\n",rt->split,rt->level,rt->type,rtlt->nn,rtgt->nn,rtlt->ntri,rtgt->ntri);
    if (rtlt->nn == 1) {
      printf("# n1 %g %g\n",n[rtlt->n[0]].c[0],n[rtlt->n[0]].c[1]);
    }
    if (rtgt->nn == 1) {
      printf("# n2 %g %g\n",n[rtgt->n[0]].c[0],n[rtgt->n[0]].c[1]);
    }
    print_rtreenode(n,tris,rtlt);
    print_rtreenode(n,tris,rtgt);
    */

    if (rtlt->level < maxlevel) {
      if (intq_pushback(&nodeq,rt->lt) != 0) {
	reportUTriError("Error: intq_pushback failed in create_utri_rtree");
	return -1;
      }
      if (intq_pushback(&nodeq,rt->gt) != 0) {
	reportUTriError("Error: intq_pushback failed in create_utri_rtree");
	return -1;
      }
    }

  }

  return 0;
}



/*
   sort nodes in order of coordinates
*/
int node_sort(node_t *a,node_t *b)
{
  if (a[0].c[0] > b[0].c[0]) return 1;
  else if (a[0].c[0] < b[0].c[0]) return -1;
  else {
    if (a[0].c[1] > b[0].c[1]) return 1;
    else if (a[0].c[1] < b[0].c[1]) return -1;
    else return 0;
  }
}

/*
   (Re)calculate a rptree for the passed table. ttype denotes the
   independent variables for the table.
   ( 0 == rho,T ; 1 == rho,E ; 2 == P,H ; 3 == P,S ; 4 == P,E )
  
   The calc_tri_rtree function assumes a single node per
   location. Thus, duplicate nodes are eliminated before the call.
*/
int calculate_rptree(eos_table_utri_t *table,int ttype)
{
  int i;
  node_t *n;
  tri_t *tris;
  rtreenodevec_t rtnv;
  int ni;
  int l;
  node_t *orign;
  int debug = 0;

  /* sanity checks */
  if (ttype < 0 || ttype > 4) {
    reportUTriError("calculate_rptree: bad table type passed");
    return -1;
  }
  if (table == NULL) {
    reportUTriError("calculate_rptree: NULL input detected");
    return -1;
  }

  /* allocate and fill node structures */
  if ((n = (node_t *)malloc(sizeof(node_t)*table->np)) == NULL) {
    reportUTriError("calculate_rptree: failed to allocate node array");
    return -1;
  }
  for (i=0;i<table->np;i++) {
    switch (ttype) {
    case 0:
      n[i].c[0] = table->modes[0]->rho[i];
      n[i].c[1] = table->modes[0]->T[i];
      break;
    case 1:
      n[i].c[0] = table->modes[0]->rho[i];
      n[i].c[1] = table->modes[0]->E[i];
      break;
    case 2:
      n[i].c[0] = table->modes[0]->p[i];
      n[i].c[1] = table->modes[0]->H[i];
      break;
    case 3:
      n[i].c[0] = table->modes[0]->p[i];
      n[i].c[1] = table->modes[0]->S[i];
      break;
    case 4:
      n[i].c[0] = table->modes[0]->p[i];
      n[i].c[1] = table->modes[0]->E[i];
      break;
    default:
      reportUTriError("calculate_rptree: bad table type in node assignment switch");
      break;
    }
    n[i].n = i;
    n[i].nn = -1;
  }

  /* sort nodes according to location */
  qsort(n,table->np,sizeof(node_t),(int (*)(const void *,const void *))node_sort);
  if (debug > 0) {
    for (i=0;i<table->np;i++) {
      printf("sn %d on %d x %g y %g\n",i,n[i].n,n[i].c[0],n[i].c[1]);
    }
  }

  /* determine mapping to other nodes and new locations */
  ni = 0;  
  n[0].nn = ni++;
  l = 0;
  for (i=1;i<table->np;i++) {
    /* check if in same location as last node */
    if (n[l].c[0] == n[i].c[0] && n[l].c[1] == n[i].c[1]) {
      n[i].nn = n[l].nn;
    }
    else {
      n[i].nn = ni++;
      l = i;
    }
  }
  if (debug > 0) {
    for (i=0;i<table->np;i++) {
      printf("sn2 %d on %d nn %d x %g y %g\n",i,n[i].n,n[i].nn,n[i].c[0],n[i].c[1]);
    }
  }

  /* unsort nodes and make new node array */
  if ((orign = (node_t *)malloc(sizeof(node_t)*table->np)) == NULL) {
    reportUTriError("calculate_rptree: failed to allocate orign array");
    return -1;
  }
  for (i=0;i<table->np;i++) {
    orign[n[i].n] = n[i];
  }
  for (i=0;i<table->np;i++) {
    n[orign[i].nn] = orign[i];
  }
  if (debug > 0) {
    for (i=0;i<table->np;i++) {
      printf("orign %d on %d nn %d x %g y %g\n",i,orign[i].n,orign[i].nn,orign[i].c[0],orign[i].c[1]);
    }
    for (i=0;i<ni;i++) {
      printf("n %d on %d nn %d x %g y %g\n",i,n[i].n,n[i].nn,n[i].c[0],n[i].c[1]);
    }
  }
  
  /* allocate and fill tri structures */
  if ((tris = (tri_t *)malloc(sizeof(tri_t)*table->nt)) == NULL) {
    reportUTriError("calculate_rptree: failed to allocate tri array");
    return -1;
  }

  for (i=0;i<table->nt;i++) {
    /* map nodes to new locations */
    tris[i].n[0] = orign[table->tris[i].p[0]].nn;
    tris[i].n[1] = orign[table->tris[i].p[1]].nn;
    tris[i].n[2] = orign[table->tris[i].p[2]].nn;
    setbb(&tris[i],n);
    if (debug > 0) printf("tri %d n %d %d %d\n",i,tris[i].n[0],tris[i].n[1],tris[i].n[2]);
  }
  if (debug > 0) {
    printf("\n\n");
    for (i=0;i<table->nt;i++) {
      print_tri(n,&tris[i]);
    }
    printf("\n\n");
  }

  /* create rtree using contracted node array */
  if (create_utri_rtree(&rtnv,n,ni,tris,table->nt) != 0) {
    reportUTriError("calculate_rptree: create_utri_rtree failed");
    return -1;
  }

  /* map tri node locations back to original values for saving */
  for (i=0;i<table->nt;i++) {
    tris[i].n[0] = table->tris[i].p[0];
    tris[i].n[1] = table->tris[i].p[1];
    tris[i].n[2] = table->tris[i].p[2];
  }

  /* save using original node array */
  if (save_utri_rtree(table,tris,table->nt,orign,table->np,&rtnv) != 0) {
    reportUTriError("calculate_rptree: save_utri_rtree failed");
    return -1;
  }

  free(orign);
  free(n);
  free(tris);
  rtreenodevec_destroy(&rtnv);

  return 0;
}
