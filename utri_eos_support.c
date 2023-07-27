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
#include <unistd.h>

#include "netcdf.h"

#include "utri_FC.h"
#include "utri_eos_support.h"
#include "ncstatus.h"

#define CHECK_NETCDF_STATUS(status,ncfun,database,fun) do { \
  if (bad_netcdf_status(status,#ncfun,database,#fun,__FILE__,__LINE__)) return status; \
} while (0)

/* type to hold information for an instance of the utri eos library */
typedef struct utri_eos_handle
{
  /* broadcast function for this instance */
  void (*bcast_fun)( void * bcast_data,
                     void * data,
                     int n,
                     int * error );

  /* state data for broadcast function
   * this must be allocated by internal routines as it is freed on cleanup */
  void *bcast_data;

  /* true (1) if this process is the I/O process */
  int master;

} utri_eos_handle_t;

/* array of handles */
int n_handles = 0;
utri_eos_handle_t * handles = NULL;

/*
 * Initialize a new handle.
 *
 */
int utri_eos_handle_init( void (*fun)( void *, void *, int, int * ),
                          void * data,
                          int master )
{
  int i,h;
  utri_eos_handle_t * tmph;

  /* find an unused handle */
  h = -1;
  for (i=0;i<n_handles;i++) {
    if (handles[i].master < 0) {
      h=i;
      break;
    }
  }

  /* reallocate handles */
  if (h < 0) {
    tmph = (utri_eos_handle_t *)realloc(handles,(n_handles+1)*sizeof(utri_eos_handle_t));
    if (tmph == NULL) {
      perror("utri_eos_handle_init: unable to reallocate handle memory");
      return -1;
    }
    h = n_handles++;
    handles = tmph;
  }

  /* save new handle */
  handles[h].bcast_fun = fun;
  handles[h].bcast_data = data;
  handles[h].master = master;

  /* return the new handle */
  return h;
}

/*
 * Close and free resources associated with a utri eos handle. Also
 * free the handles array if none in use.
 *
 */
void utri_eos_handle_free( const int handle )
{
  int i;

  /* skip invalid handles */
  if (handle < 0 || handle >= n_handles) return;

  /* clean out this handle */
  handles[handle].bcast_fun = NULL;
  if (handles[handle].bcast_data) free(handles[handle].bcast_data);
  handles[handle].bcast_data = NULL;
  handles[handle].master = -1;

  /* free the handles array once all are empty. If this is not the
   * case, meaning master != -1 for some entry, then skip clean up.
   */
  for (i=0;i<n_handles;i++) if (handles[i].master != -1) return;

  /* okay to clean up */
  if (handles) free(handles);
  handles = NULL;
}

/*
 * Generic data broadcast function.
 *
 * inputs:
 *   handle -- index to the utri eos handles
 *   data   -- array of data (significant only on master node)
 *   n      -- size of the data to be broadcast
 *
 * outputs:
 *   data  -- broadcast array of data
 *   error -- error status (non-zero for a broadcast failure)
 *
 */
void broadcast_utri_data( int handle,
                          void * data,
                          int n,
                          int * error )
{
  if (handles[handle].bcast_fun)
    handles[handle].bcast_fun(handles[handle].bcast_data,data,n,error);
  else *error = 0;
}

/* query the master flag for this handle */
int is_master_node( int handle )
{
  return handles[handle].master;
}

/* current error string */
char * utri_error_string = 0;

/* return the error string and clear the internal buffer */
char * utri_error_string_get()
{
  char * estring = utri_error_string;
  utri_error_string = 0;
  return estring;
}

/* append to the current error string */
void reportUTriError( const char * message )
{
  char * newmessage;
  size_t len;

  /* message length + null + newline + current length */
  len = strlen(message)+2;
  if (utri_error_string != 0) len += strlen(utri_error_string);

  /* allocate the new message string */
  newmessage = (char *)realloc(utri_error_string,len*sizeof(char));
  if (newmessage == NULL) {
    perror("Error: unable to allocate memory in reportUTriError");
    return;
  }
  else {
    if (utri_error_string == 0) newmessage[0] = '\0';
    utri_error_string = newmessage;
  }
  strcat(utri_error_string, message);
  strcat(utri_error_string, "\n");
}

typedef struct point_array_metadata {
  char *name;
  nc_type type;
  int varid;
  void *data;
} pam_t;

enum array_index_names {
  PHASE,
  R,
  T,
  P,
  E,
  S,
  F,
  DPDR,
  DPDT,
  DEDR,
  DEDT,
  CS,
  DCSDR,
  VISC,
  THCON,
  NUM_ARRAY_INDEX_NAMES
};

void load_array_metadata(pam_t *md,eos_data_points_t *edp)
{
  md[PHASE] = (pam_t) {"PHASE", NC_INT,    0, edp->phase};
  md[R]     = (pam_t) {"R",     NC_DOUBLE, 0, edp->rho  };
  md[T]     = (pam_t) {"T",     NC_DOUBLE, 0, edp->T    };
  md[P]     = (pam_t) {"P",     NC_DOUBLE, 0, edp->p    };
  md[E]     = (pam_t) {"E",     NC_DOUBLE, 0, edp->E    };
  md[S]     = (pam_t) {"S",     NC_DOUBLE, 0, edp->S    };
  md[F]     = (pam_t) {"F",     NC_DOUBLE, 0, edp->F    };
  md[DPDR]  = (pam_t) {"DPDR",  NC_DOUBLE, 0, edp->dpdr };
  md[DPDT]  = (pam_t) {"DPDT",  NC_DOUBLE, 0, edp->dpdT };
  md[DEDR]  = (pam_t) {"DEDR",  NC_DOUBLE, 0, edp->dEdr };
  md[DEDT]  = (pam_t) {"DEDT",  NC_DOUBLE, 0, edp->dEdT };
  md[CS]    = (pam_t) {"CS",    NC_DOUBLE, 0, edp->cs   };
  md[DCSDR] = (pam_t) {"DCSDR", NC_DOUBLE, 0, edp->dcsdr};
  md[VISC]  = (pam_t) {"VISC",  NC_DOUBLE, 0, edp->visc };
  md[THCON] = (pam_t) {"THCON", NC_DOUBLE, 0, edp->thcon};
}

#define MAX_NAME_LENGTH 20

int eostableutriwritenc(const char *nc_name,eos_table_utri_t *t)
{
  int status;
  int ncid;
  int points_dim;
  int modes_dim;
  int i,j;
  pam_t pa_md[NUM_ARRAY_INDEX_NAMES];
  eos_data_points_t *edp;
  int *intnodedata;
  int nind;
  int ni_dim,nd_dim;
  int ni_id,nd_id;
  int fail;
  char buff[MAX_NAME_LENGTH];

  /*
     create file
  */
  status = nc_create(nc_name,NC_CLOBBER,&ncid);
  CHECK_NETCDF_STATUS(status,nc_create,nc_name,eostableutriwritenc);

  /*
     write defaults
  */
  status = nc_def_dim(ncid,"num_defaults",5,&nd_dim);
  CHECK_NETCDF_STATUS(status,nc_def_dim,nc_name,eostableutriwritenc);

  status = nc_def_var(ncid,"defaults",NC_DOUBLE,1,&nd_dim,&nd_id);
  CHECK_NETCDF_STATUS(status,nc_def_var,nc_name,eostableutriwritenc);

  status = nc_enddef(ncid);
  CHECK_NETCDF_STATUS(status,nc_enddef,nc_name,eostableutriwritenc);

  status = nc_put_var_double(ncid,nd_id,t->defaults);
  CHECK_NETCDF_STATUS(status,nc_put_var_double,nc_name,eostableutriwritenc);

  /*
     write flags
  */
  status = nc_redef(ncid);
  CHECK_NETCDF_STATUS(status,nc_redef,nc_name,eostableutriwritenc);

  status = nc_def_dim(ncid,"num_flags",2,&ni_dim);
  CHECK_NETCDF_STATUS(status,nc_def_dim,nc_name,eostableutriwritenc);

  status = nc_def_var(ncid,"flags",NC_INT,1,&ni_dim,&ni_id);
  CHECK_NETCDF_STATUS(status,nc_def_var,nc_name,eostableutriwritenc);

  status = nc_enddef(ncid);
  CHECK_NETCDF_STATUS(status,nc_enddef,nc_name,eostableutriwritenc);

  status = nc_put_var_int(ncid,ni_id,t->flags);
  CHECK_NETCDF_STATUS(status,nc_put_var_int,nc_name,eostableutriwritenc);

  /*
     write boundary information
  */
  if ((intnodedata = (int *)malloc(sizeof(int)*6)) == NULL) {
    printf("failed to allocate memory in eostableutriwritenc\n");
    return 1;
  }
  intnodedata[0] = t->nbp[0];
  intnodedata[1] = t->nbp[1];
  intnodedata[2] = t->nbt[0];
  intnodedata[3] = t->nbt[1];
  intnodedata[4] = t->nbi[0];
  intnodedata[5] = t->nbi[1];

  status = nc_redef(ncid);
  CHECK_NETCDF_STATUS(status,nc_redef,nc_name,eostableutriwritenc);

  status = nc_def_dim(ncid,"num_bounds",6,&ni_dim);
  CHECK_NETCDF_STATUS(status,nc_def_dim,nc_name,eostableutriwritenc);

  status = nc_def_var(ncid,"bounds",NC_INT,1,&ni_dim,&ni_id);
  CHECK_NETCDF_STATUS(status,nc_def_var,nc_name,eostableutriwritenc);

  status = nc_enddef(ncid);
  CHECK_NETCDF_STATUS(status,nc_enddef,nc_name,eostableutriwritenc);

  status = nc_put_var_int(ncid,ni_id,intnodedata);
  CHECK_NETCDF_STATUS(status,nc_put_var_int,nc_name,eostableutriwritenc);

  free(intnodedata);

  /*
     write modes
  */
  if (t->nmodes < 1) {
    printf("eostableutriwritenc: modes length less than 1 in table\n");
    return 1;
  }

  for (i=0;i<t->nmodes;i++) {
    edp = (eos_data_points_t *) t->modes[i];
    if (t->np != edp->n) {
      printf("eostableutriwritenc: mode %d points length inconsistent in utri %d %d\n",i,t->np,edp->n);
      return 1;
    }
  }

  status = nc_redef(ncid);
  CHECK_NETCDF_STATUS(status,nc_redef,nc_name,eostableutriwritenc);

  status = nc_def_dim(ncid,"num_points",t->np,&points_dim);
  CHECK_NETCDF_STATUS(status,nc_def_dim,nc_name,eostableutriwritenc);

  status = nc_def_dim(ncid,"num_modes",t->nmodes,&modes_dim);
  CHECK_NETCDF_STATUS(status,nc_def_dim,nc_name,eostableutriwritenc);

  status = nc_enddef(ncid);
  CHECK_NETCDF_STATUS(status,nc_enddef,nc_name,eostableutriwritenc);

  for (j=0;j<t->nmodes;j++) {

    load_array_metadata(pa_md,(eos_data_points_t *) t->modes[j]);

    status = nc_redef(ncid);
    CHECK_NETCDF_STATUS(status,nc_redef,nc_name,eostableutriwritenc);

    for (i=0;i<NUM_ARRAY_INDEX_NAMES;i++) {
      if (snprintf(buff,MAX_NAME_LENGTH,
                   "%s%d",pa_md[i].name,j) >= MAX_NAME_LENGTH) {
        printf("eostableutriwritenc: string of length %d insufficient to hold mode variable name %s%d\n",MAX_NAME_LENGTH,pa_md[i].name,j);
        return 1;
      }
      status = nc_def_var(ncid,buff,pa_md[i].type,1,&points_dim,&pa_md[i].varid);
      CHECK_NETCDF_STATUS(status,nc_def_var,nc_name,eostableutriwritenc);
    }

    status = nc_enddef(ncid);
    CHECK_NETCDF_STATUS(status,nc_enddef,nc_name,eostableutriwritenc);

    for (i=0;i<NUM_ARRAY_INDEX_NAMES;i++) {
      if (pa_md[i].type == NC_INT) {
        status = nc_put_var_int(ncid,pa_md[i].varid,(int *)pa_md[i].data);
        CHECK_NETCDF_STATUS(status,nc_put_var_int,nc_name,eostableutriwritenc);
      }
      else if (pa_md[i].type == NC_DOUBLE) {
        status = nc_put_var_double(ncid,pa_md[i].varid,(double *)pa_md[i].data);
        CHECK_NETCDF_STATUS(status,nc_put_var_double,nc_name,eostableutriwritenc);
      }
      else {
        printf("unsupported variable type in eostableutriwritenc\n");
        return 1;
      }
    }
  }

  /*
     write tri nodes
  */
  nind = t->nt*3;
  fail = 0;
  if ((intnodedata = (int *)malloc(sizeof(int)*nind)) == NULL) fail = 1;
  if (fail) {
    printf("failed to allocate memory in eostableutriwritenc\n");
    return 1;
  }
  for (i=0;i<t->nt;i++) {
    intnodedata[3*i+0] = t->tris[i].p[0];
    intnodedata[3*i+1] = t->tris[i].p[1];
    intnodedata[3*i+2] = t->tris[i].p[2];
  }

  status = nc_redef(ncid);
  CHECK_NETCDF_STATUS(status,nc_redef,nc_name,eostableutriwritenc);

  status = nc_def_dim(ncid,"num_tri_ints",nind,&ni_dim);
  CHECK_NETCDF_STATUS(status,nc_def_dim,nc_name,eostableutriwritenc);

  status = nc_def_var(ncid,"tri_ints",NC_INT,1,&ni_dim,&ni_id);
  CHECK_NETCDF_STATUS(status,nc_def_var,nc_name,eostableutriwritenc);

  status = nc_enddef(ncid);
  CHECK_NETCDF_STATUS(status,nc_enddef,nc_name,eostableutriwritenc);

  status = nc_put_var_int(ncid,ni_id,intnodedata);
  CHECK_NETCDF_STATUS(status,nc_put_var_int,nc_name,eostableutriwritenc);

  free(intnodedata);

  /*
     close file
  */
  status = nc_close(ncid);
  CHECK_NETCDF_STATUS(status,nc_close,nc_name,eostableutriwritenc);

  return NC_NOERR;
}

int eostableutrireadnc(const char *nc_name,eos_table_utri_t *t, int handle)
{
  int i,j;
  int status;
  int ncid;
  int points_dim;
  int modes_dim;
  nc_type point_type;
  int point_dim;
  int point_ndim;
  size_t length;
  size_t totmodes;
  eos_data_points_t *edp=NULL;
  pam_t pa_md[NUM_ARRAY_INDEX_NAMES];
  int ni_dim;
  int ni_id;
  size_t nind;
  int fail;
  int *intnodedata = NULL;
  double *doublenodedata = NULL;
  int defaults_dim;
  int defaults_id;
  char buff[MAX_NAME_LENGTH];
  char mesg[80];
  int master = is_master_node(handle);
  int error = 0;
  int flags_dim;
  int flags_id;

  /*
     open file
  */
  if (master) {
    status = nc_open(nc_name,NC_NOWRITE,&ncid);
    CHECK_NETCDF_STATUS(status,nc_open,nc_name,eostableutrireadnc);

    /*
       load flags
    */
    status = nc_inq_dimid(ncid,"num_flags",&flags_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);

    status = nc_inq_dimlen(ncid,flags_dim,&length);
    CHECK_NETCDF_STATUS(status,nc_inq_dimlen,nc_name,eostableutrireadnc);
    if (length != 2) {
      snprintf(mesg,80,"invalid flags array length (%d) reading %s in eostableutrireadnc",
         (int)length,nc_name);
      reportUTriError(mesg);
      return 1;
    }

    status = nc_inq_varid(ncid,"flags",&flags_id);
    CHECK_NETCDF_STATUS(status,nc_inq_varid,nc_name,eostableutrireadnc);
    status = nc_inq_vartype(ncid,flags_id,&point_type);
    CHECK_NETCDF_STATUS(status,nc_inq_vartype,nc_name,eostableutrireadnc);
    if (point_type != NC_INT) {
      snprintf(mesg,80,"bad point type for flags reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_varndims(ncid,flags_id,&point_ndim);
    CHECK_NETCDF_STATUS(status,nc_inq_varndims,nc_name,eostableutrireadnc);
    if (point_ndim != 1) {
      snprintf(mesg,80,"flags has wrong number of dimensions reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_vardimid(ncid,flags_id,&point_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);
    if (point_dim != flags_dim) {
      snprintf(mesg,80,"flags has incorrect dimension reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_get_var_int(ncid,flags_id,t->flags);
    CHECK_NETCDF_STATUS(status,nc_get_var_int,nc_name,eostableutrireadnc);

    /*
       load defaults
    */
    status = nc_inq_dimid(ncid,"num_defaults",&defaults_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);

    status = nc_inq_dimlen(ncid,defaults_dim,&length);
    CHECK_NETCDF_STATUS(status,nc_inq_dimlen,nc_name,eostableutrireadnc);
    if (length != 5) {
      snprintf(mesg,80,"invalid defaults array length (%d) reading %s in eostableutrireadnc",
         (int)length,nc_name);
      reportUTriError(mesg);
      return 1;
    }
  }
  broadcast_utri_data(handle, t->flags, sizeof(int)*2, &error);
  broadcast_utri_data(handle, &defaults_dim, sizeof(int), &error);
  broadcast_utri_data(handle, &length, sizeof(size_t), &error);

  if ((doublenodedata = (double *)malloc(sizeof(double)*length)) == NULL) {
    reportUTriError("failed to allocate memory in eostableutrireadnc");
    return 1;
  }

  if (master) {
    status = nc_inq_varid(ncid,"defaults",&defaults_id);
    CHECK_NETCDF_STATUS(status,nc_inq_varid,nc_name,eostableutrireadnc);
    status = nc_inq_vartype(ncid,defaults_id,&point_type);
    CHECK_NETCDF_STATUS(status,nc_inq_vartype,nc_name,eostableutrireadnc);
    if (point_type != NC_DOUBLE) {
      snprintf(mesg,80,"bad point type for defaults reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_varndims(ncid,defaults_id,&point_ndim);
    CHECK_NETCDF_STATUS(status,nc_inq_varndims,nc_name,eostableutrireadnc);
    if (point_ndim != 1) {
      snprintf(mesg,80,"defaults has wrong number of dimensions reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_vardimid(ncid,defaults_id,&point_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);
    if (point_dim != defaults_dim) {
      snprintf(mesg,80,"defaults has incorrect dimension reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_get_var_double(ncid,defaults_id,doublenodedata);
    CHECK_NETCDF_STATUS(status,nc_get_var_double,nc_name,eostableutrireadnc);
  }
  broadcast_utri_data(handle, &point_type, sizeof(int), &error);
  broadcast_utri_data(handle, &point_ndim, sizeof(int), &error);
  broadcast_utri_data(handle, &point_dim, sizeof(int), &error);
  broadcast_utri_data(handle, doublenodedata, sizeof(double)*length, &error);

  t->defaults = doublenodedata;

  /*
     load boundary information
  */
  if (master) {
    status = nc_inq_dimid(ncid,"num_bounds",&ni_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);

    status = nc_inq_dimlen(ncid,ni_dim,&length);
    CHECK_NETCDF_STATUS(status,nc_inq_dimlen,nc_name,eostableutrireadnc);
    if (length != 6) {
      snprintf(mesg,80,"invalid bounds array length (%d) reading %s in eostableutrireadnc",
         (int)length,nc_name);
      reportUTriError(mesg);
      return 1;
    }
  }
  broadcast_utri_data(handle, &ni_dim, sizeof(int), &error);
  broadcast_utri_data(handle, &length, sizeof(size_t), &error);

  if ((intnodedata = (int *)malloc(sizeof(int)*length)) == NULL) {
    reportUTriError("failed to allocate memory in eostableutrireadnc");
    return 1;
  }

  if (master) {
    status = nc_inq_varid(ncid,"bounds",&ni_id);
    CHECK_NETCDF_STATUS(status,nc_inq_varid,nc_name,eostableutrireadnc);
    status = nc_inq_vartype(ncid,ni_id,&point_type);
    CHECK_NETCDF_STATUS(status,nc_inq_vartype,nc_name,eostableutrireadnc);
    if (point_type != NC_INT) {
      snprintf(mesg,80,"bad point type for bounds reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_varndims(ncid,ni_id,&point_ndim);
    CHECK_NETCDF_STATUS(status,nc_inq_varndims,nc_name,eostableutrireadnc);
    if (point_ndim != 1) {
      snprintf(mesg,80,"bounds has wrong number of dimensions reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_vardimid(ncid,ni_id,&point_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);
    if (point_dim != ni_dim) {
      snprintf(mesg,80,"bounds has incorrect dimension reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_get_var_int(ncid,ni_id,intnodedata);
    CHECK_NETCDF_STATUS(status,nc_get_var_int,nc_name,eostableutrireadnc);
  }
  broadcast_utri_data(handle, &ni_id, sizeof(int), &error);
  broadcast_utri_data(handle, &point_type, sizeof(int), &error);
  broadcast_utri_data(handle, &point_ndim, sizeof(int), &error);
  broadcast_utri_data(handle, &point_dim, sizeof(int), &error);
  broadcast_utri_data(handle, intnodedata, sizeof(int)*length, &error);

  t->nbp[0] = intnodedata[0];
  t->nbp[1] = intnodedata[1];
  t->nbt[0] = intnodedata[2];
  t->nbt[1] = intnodedata[3];
  t->nbi[0] = intnodedata[4];
  t->nbi[1] = intnodedata[5];

  free(intnodedata);

  /*
     load modes
  */
  if (master) {
    status = nc_inq_dimid(ncid,"num_points",&points_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);

    status = nc_inq_dimlen(ncid,points_dim,&length);
    CHECK_NETCDF_STATUS(status,nc_inq_dimlen,nc_name,eostableutrireadnc);

    if (length < 1) {
      snprintf(mesg,80,"invalid point array length (%d) reading %s in eostableutrireadnc",
         (int)length,nc_name);
      reportUTriError(mesg);
      return 1;
    }
  }
  broadcast_utri_data(handle, &point_dim, sizeof(int), &error);
  broadcast_utri_data(handle, &length, sizeof(size_t), &error);

  t->np = length;

  if (master) {
    status = nc_inq_dimid(ncid,"num_modes",&modes_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);

    status = nc_inq_dimlen(ncid,modes_dim,&totmodes);
    CHECK_NETCDF_STATUS(status,nc_inq_dimlen,nc_name,eostableutrireadnc);

    if (totmodes < 1) {
      snprintf(mesg,80,"invalid number of modes (%d) reading %s in eostableutrireadnc",
         (int)totmodes,nc_name);
      reportUTriError(mesg);
      return 1;
    }
  }
  broadcast_utri_data(handle, &modes_dim, sizeof(int), &error);
  broadcast_utri_data(handle, &totmodes, sizeof(size_t), &error);

  t->nmodes = totmodes;
  if ((t->modes = (eos_data_points_t **)malloc(sizeof(eos_data_points_t *)*totmodes)) == NULL) {
    reportUTriError("failed to allocate mode array in eostableutrireadnc");
    return 1;
  }

  for (j=0;j<t->nmodes;j++) {

    edp = NULL;
    if ((edp = (eos_data_points_t *) eospointsrealloc((void *)edp,length)) == NULL) {
      snprintf(mesg,80,"failed to allocate memory for mode %d edp in eostableutrireadnc",j);
      reportUTriError(mesg);
      return 1;
    }

    load_array_metadata(pa_md,edp);

    for (i=0;i<NUM_ARRAY_INDEX_NAMES;i++) {
      if (master) {
        if (snprintf(buff,MAX_NAME_LENGTH,
         "%s%d",pa_md[i].name,j) >= MAX_NAME_LENGTH) {
          snprintf(mesg,80,"eostableutrireadnc: string of length %d insufficient to hold mode variable name %s%d",MAX_NAME_LENGTH,pa_md[i].name,j);
          reportUTriError(mesg);
          return 1;
        }
        status = nc_inq_varid(ncid,buff,&pa_md[i].varid);
        CHECK_NETCDF_STATUS(status,nc_inq_varid,nc_name,eostableutrireadnc);
        status = nc_inq_vartype(ncid,pa_md[i].varid,&point_type);
        CHECK_NETCDF_STATUS(status,nc_inq_vartype,nc_name,eostableutrireadnc);
        if (point_type != pa_md[i].type) {
          snprintf(mesg,80,"bad point type for %s reading %s in eostableutrireadnc",
            buff,nc_name);
          reportUTriError(mesg);
          return 1;
        }
        status = nc_inq_varndims(ncid,pa_md[i].varid,&point_ndim);
        CHECK_NETCDF_STATUS(status,nc_inq_varndims,nc_name,eostableutrireadnc);
        if (point_ndim != 1) {
          snprintf(mesg,80,"%s has wrong number of dimensions reading %s in eostableutrireadnc",
            buff,nc_name);
          reportUTriError(mesg);
          return 1;
        }
        status = nc_inq_vardimid(ncid,pa_md[i].varid,&point_dim);
        CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);
        if (point_dim != points_dim) {
          snprintf(mesg,80,"%s has incorrect dimension reading %s in eostableutrireadnc",
            buff,nc_name);
          reportUTriError(mesg);
          return 1;
        }
        if (pa_md[i].type == NC_INT) {
          status = nc_get_var_int(ncid,pa_md[i].varid,(int *)pa_md[i].data);
          CHECK_NETCDF_STATUS(status,nc_get_var_int,nc_name,eostableutrireadnc);
        }
        else if (pa_md[i].type == NC_DOUBLE) {
          status = nc_get_var_double(ncid,pa_md[i].varid,(double *)pa_md[i].data);
          CHECK_NETCDF_STATUS(status,nc_get_var_double,nc_name,eostableutrireadnc);
        }
        else {
          reportUTriError("unsupported variable type in eostableutrireadnc");
          return 1;
        }
      }
      broadcast_utri_data(handle, &pa_md[i].varid, sizeof(int), &error);
      broadcast_utri_data(handle, &point_type, sizeof(int), &error);
      broadcast_utri_data(handle, &point_ndim, sizeof(int), &error);
      broadcast_utri_data(handle, &point_dim, sizeof(int), &error);
      if (pa_md[i].type == NC_INT)
        broadcast_utri_data(handle, pa_md[i].data, sizeof(int)*edp->n, &error);
      else if (pa_md[i].type == NC_DOUBLE)
        broadcast_utri_data(handle, pa_md[i].data, sizeof(double)*edp->n, &error);
    }

    t->modes[j] = edp;
  }

  /*
     load tris
  */
  if (master) {
    status = nc_inq_dimid(ncid,"num_tri_ints",&ni_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);

    status = nc_inq_dimlen(ncid,ni_dim,&nind);
    CHECK_NETCDF_STATUS(status,nc_inq_dimlen,nc_name,eostableutrireadnc);

    if (nind < 1 || nind%3 != 0) {
      snprintf(mesg,80,"invalid tri ints array length (%d) reading %s in eostableutrireadnc",
         (int)nind,nc_name);
      reportUTriError(mesg);
      return 1;
    }
  }
  broadcast_utri_data(handle, &ni_dim, sizeof(int), &error);
  broadcast_utri_data(handle, &nind, sizeof(size_t), &error);

  t->nt = nind/3;

  fail = 0;
  if ((intnodedata = (int *)malloc(sizeof(int)*nind)) == NULL) fail = 1;
  if (fail) {
    reportUTriError("failed to allocate memory in eostableutrireadnc\n");
    return 1;
  }

  if (master) {
    status = nc_inq_varid(ncid,"tri_ints",&ni_id);
    CHECK_NETCDF_STATUS(status,nc_inq_varid,nc_name,eostableutrireadnc);
    status = nc_inq_vartype(ncid,ni_id,&point_type);
    CHECK_NETCDF_STATUS(status,nc_inq_vartype,nc_name,eostableutrireadnc);
    if (point_type != NC_INT) {
      snprintf(mesg,80,"bad point type for tri_ints reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_varndims(ncid,ni_id,&point_ndim);
    CHECK_NETCDF_STATUS(status,nc_inq_varndims,nc_name,eostableutrireadnc);
    if (point_ndim != 1) {
      snprintf(mesg,80,"tri_ints has wrong number of dimensions reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_inq_vardimid(ncid,ni_id,&point_dim);
    CHECK_NETCDF_STATUS(status,nc_inq_dimid,nc_name,eostableutrireadnc);
    if (point_dim != ni_dim) {
      snprintf(mesg,80,"tri_ints has incorrect dimension reading %s in eostableutrireadnc",nc_name);
      reportUTriError(mesg);
      return 1;
    }
    status = nc_get_var_int(ncid,ni_id,intnodedata);
    CHECK_NETCDF_STATUS(status,nc_get_var_int,nc_name,eostableutrireadnc);
  }
  broadcast_utri_data(handle, &ni_id, sizeof(int), &error);
  broadcast_utri_data(handle, &point_type, sizeof(int), &error);
  broadcast_utri_data(handle, &point_ndim, sizeof(int), &error);
  broadcast_utri_data(handle, &point_dim, sizeof(int), &error);
  broadcast_utri_data(handle, intnodedata, sizeof(int)*nind, &error);

  if ((t->tris = (eos_table_tri_t *)malloc(sizeof(eos_table_tri_t)*t->nt)) == NULL) {
    reportUTriError("failed to allocate memory in eostableutrireadnc");
    return 1;
  }

  for (i=0;i<t->nt;i++) {
    t->tris[i].p[0] = intnodedata[3*i+0];
    t->tris[i].p[1] = intnodedata[3*i+1];
    t->tris[i].p[2] = intnodedata[3*i+2];
    t->tris[i].a[0] = 0.;
    t->tris[i].a[1] = 0.;
    t->tris[i].a[2] = 0.;
    t->tris[i].a[3] = 0.;
    t->tris[i].a[4] = 0.;
    t->tris[i].a[5] = 0.;
  }

  free(intnodedata);

  t->rptree.splits = NULL;
  t->rptree.trilookup = NULL;
  t->rptree.trilist = NULL;

  /*
     close file
  */
  if (master) {
    status = nc_close(ncid);
    CHECK_NETCDF_STATUS(status,nc_close,nc_name,eostableutrireadnc);
  }

  return NC_NOERR;
}

/*
   linear interpolation on a triangle
  
   inputs are:
    n -- number of points
    tris -- array of triangles
    itri -- array of indices to interpolation triangle
    points -- array of data points
    b -- array of barycentric coordinates
   output is:
    f -- array of interpolated values
*/
void eostableutric0interp(int n,eos_table_tri_t *tris,int *itri,double *points,double *b,double *f)
{
  int i;

  for (i=0;i<n;i++) {
    f[i] =
      b[3*i+0]*points[tris[-itri[i]].p[0]] +
      b[3*i+1]*points[tris[-itri[i]].p[1]] +
      b[3*i+2]*points[tris[-itri[i]].p[2]];
  }

}

void bc2dtcoord2(double *a,double px,double py,double *l)
{
  l[0] = a[1]*(a[5]-py)-(a[2]-px)*a[4];
  l[1] = a[3]*(a[2]-px)-(a[5]-py)*a[0];
  l[2] = 1.-l[0]-l[1];
  if (fabs(l[0]) < 2.e-14) l[0] = 0.;
  if (fabs(l[1]) < 2.e-14) l[1] = 0.;
  if (fabs(l[2]) < 2.e-14) l[2] = 0.;
}

/*
 * Vectorized entry point for evaluating functions on the unstructured triangle tables
 *
 * t       -- input table
 * stride  --  size of the dataset
 * n       -- number of entries to evaluate
 * x, y    -- input independent variable arrays
 * f       -- arrays in which to place the dependent variable output
 * type    -- denotes the set of independent variables
 *            { 1 == (rho,T)  2 == (rho,E)  3 == (P,H)  4 == (P,S)  5 == (P,E) }
 * clip    -- length 4 array on output incremented by the number of variables clipped
 * scratch -- stride*5 length array of scratch space to aid vectorization
 * int_scratch -- stride length array of integer scratch space
 *
 */
int eostableutrievaluatec0(eos_table_utri_t *t,int stride,int n,double *x,double *y,eos_data_points_t *f,int type,int *clip,double *scratch,int *int_scratch)
{
  int *nt;
  double *b;
  int i,j;
  int done;
  double *xyc;
  int xyvar;
  int offset;

  /* tri node pointer */
  nt = &int_scratch[0];
  /* barycentric coordinates pointer (3*stride length) */
  b = &scratch[0];
  /* clipped values of independent variables (2*stride length) */
  xyc = &scratch[3*stride];

  /* initialize arrays */
  for (j=0;j<n;j++) {

    /* set independent variables */
    xyc[2*j+0] = x[j];
    xyc[2*j+1] = y[j];

    /* clip out of table input values */
    if (xyc[2*j] < t->rptree.xbounds[0]) {
      xyc[2*j] = t->rptree.xbounds[0];
      clip[0]++;
    }
    else if (xyc[2*j] > t->rptree.xbounds[1]) {
      xyc[2*j] = t->rptree.xbounds[1];
      clip[1]++;
    }
    if (xyc[2*j+1] < t->rptree.ybounds[0]) {
      xyc[2*j+1] = t->rptree.ybounds[0];
      clip[2]++;
    }
    else if (xyc[2*j+1] > t->rptree.ybounds[1]) {
      xyc[2*j+1] = t->rptree.ybounds[1];
      clip[3]++;
    }

    /* initial triangle index location */
    nt[j] = 0;
  }

  /* lookup triangle list */
  xyvar = 0;
  offset = 1;
  for (i=0;i<t->rptree.slevels;i++) {
    for (j=0;j<n;j++) {
      if (xyc[2*j+xyvar] < t->rptree.splits[nt[j]]) {
        nt[j] = 2*nt[j]+1;
      }
      else {
        nt[j] = 2*nt[j]+2;
      }
    }
    xyvar = (xyvar+1)%2;
    offset *= 2;
  }

  /* adjust nt to index the trilookup array */
  for (j=0;j<n;j++) {
    nt[j] = 2*(nt[j] - offset + 1);
  }

  /* find triangle for this point */
  done = 0;
  offset = 0;
  while (done == 0) {
    done = 1;
    for (j=0;j<n;j++) {
      if (nt[j] >= 0) {
        bc2dtcoord2(t->tris[t->rptree.trilist[t->rptree.trilookup[nt[j]+1]+offset]].a,xyc[2*j],xyc[2*j+1],&b[3*j]);
        if ((b[3*j+0] >= 0. && b[3*j+0] <= 1. &&
             b[3*j+1] >= 0. && b[3*j+1] <= 1. &&
             b[3*j+2] >= 0. && b[3*j+2] <= 1.) ||
             offset+1 >= t->rptree.trilookup[nt[j]]) {
          nt[j] = -t->rptree.trilist[t->rptree.trilookup[nt[j]+1]+offset];
        }
        else done = 0;
      }
    }
    offset++;
  }

  /* interpolate */
  if (type == 1 || type == 2) {
    if (type == 1) {
      eostableutric0interp(n,t->tris,nt,t->modes[0]->E,b,f->E);
    }
    else if (type == 2) {
      eostableutric0interp(n,t->tris,nt,t->modes[0]->T,b,f->T);
    }
    /*
     * init routine has to protect from this code path
    else {
      reportUTriError("bad type in eostableutrievaluatec0");
      return -1;
    }
    */
    eostableutric0interp(n,t->tris,nt,t->modes[0]->p,b,f->p);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->cs,b,f->cs);
    for (j=0;j<n;j++) {
      if (f->cs[j] < t->defaults[4]) f->cs[j] = t->defaults[4];
    }
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dpdr,b,f->dpdr);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dpdT,b,f->dpdT);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dEdr,b,f->dEdr);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dEdT,b,f->dEdT);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dcsdr,b,f->dcsdr);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->visc,b,f->visc);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->thcon,b,f->thcon);
    for (j=0;j<n;j++) {
      f->phase[j] = t->modes[0]->phase[t->tris[-nt[j]].p[0]];
    }
  }
  else { /* type == 3-5 */
    if (type == 3) {
      eostableutric0interp(n,t->tris,nt,t->modes[0]->E,b,f->E);
      eostableutric0interp(n,t->tris,nt,t->modes[0]->S,b,f->S);
    }
    else if (type == 4) {
      eostableutric0interp(n,t->tris,nt,t->modes[0]->E,b,f->E);
      eostableutric0interp(n,t->tris,nt,t->modes[0]->dcsdr,b,f->dcsdr); /* really H */
    }
    else if (type == 5) {
      eostableutric0interp(n,t->tris,nt,t->modes[0]->dcsdr,b,f->dcsdr); /* really H */
      eostableutric0interp(n,t->tris,nt,t->modes[0]->S,b,f->S);
    }
    /*
     * init routine has to protect from this code path
    else {
      reportUTriError("bad type in eostableutrievaluatec0");
      return -1;
    }
    */
    eostableutric0interp(n,t->tris,nt,t->modes[0]->rho,b,f->rho);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->T,b,f->T);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dEdr,b,f->dEdr); /* really G */
    eostableutric0interp(n,t->tris,nt,t->modes[0]->F,b,f->F);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->cs,b,f->cs);
    for (j=0;j<n;j++) {
      if (f->cs[j] < t->defaults[4]) f->cs[j] = t->defaults[4];
    }
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dpdr,b,f->dpdr);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dpdT,b,f->dpdT);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->dEdT,b,f->dEdT);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->visc,b,f->visc);
    eostableutric0interp(n,t->tris,nt,t->modes[0]->thcon,b,f->thcon);
    for (j=0;j<n;j++) {
      f->phase[j] = t->modes[0]->phase[t->tris[-nt[j]].p[0]];
    }
  }

  return 0;
}

void eos_table_utri_delete(eos_table_utri_t *t)
{
  int i;

  if (t->rptree.splits) free(t->rptree.splits);
  if (t->rptree.trilookup) free(t->rptree.trilookup);
  if (t->rptree.trilist) free(t->rptree.trilist);
  if (t->modes) {
    for (i=0;i<t->nmodes;i++) {
      eospointsfree(t->modes[i]);
      if (t->modes[i]) free(t->modes[i]);
    }
    free(t->modes);
  }
  if (t->tris) free(t->tris);
  if (t->defaults) free(t->defaults);
}

/*
 * Construct the boundary triangles for the table
 *
 * t     -- input table
 * ttype -- type of table (2 == RE)
 *                        (5 == PE)
 *
 * For the upper and lower boundary the algorithm searches for the
 * extreme point, then sets the new table boundary slightly outside
 * this value. The boundary points are then moved to the new boundary
 * value and triangle nodes adjusted accordingly. Only mode 0 is
 * modified in the table.
 *
 */
int construct_boundary_tris(eos_table_utri_t *t,int ttype)
{
  double dir;
  int b;
  int i;
  double maxb;
  int ii,iin;
  int tritype;
  int tmp;
  int tri;

  /*
     only handle energy tables for now
  */
  if (ttype != 2 && ttype != 5) {
    reportUTriError("Error in construct_boundary_tris: only table type 2 supported");
    return -1;
  }

  /*
     check that there is a boundary on this table
  */
  if (t->nbp[0] < 1 || t->nbp[1] < 1) {
    reportUTriError("Error in construct_boundary_tris: table missing boundary nodes");
    return -1;
  }

  /*
     construct the boundary triangles for each boundary
  */
  for (b=0;b<2;b++) {
    if (b == 0) dir = -1.;
    else dir = 1.;

    /* find max value */
    maxb = -9.e99*dir;
    for (i=0;i<t->nbp[b];i++) {
      ii = t->nbi[b]+i;
      if (ttype == 2 || ttype == 5) {
        if (t->modes[0]->E[ii]*dir > maxb*dir) maxb = t->modes[0]->E[ii];
      }
    }

    /* change boundary points to have more than max value */
    for (i=0;i<t->nbp[b];i++) {
      ii = t->nbi[b]+i;
      if (ttype == 2 || ttype == 5) t->modes[0]->E[ii] = maxb+dir*0.01*(1.+fabs(maxb));
    }

    /* adjust triangles */
    tri = t->nbt[b];
    for (i=0;i<t->nbp[b]-1;i++) {
      ii = t->nbi[b]+i;
      iin = ii + t->nbp[b];

      /* only make tris with the same phase */
      if (t->modes[0]->phase[ii] != t->modes[0]->phase[ii+1]) continue;

      /* choose the nodes to link in the tri */
      /*if (ttype == 2 || ttype == 5) { */
        if (t->modes[0]->E[iin]*dir > t->modes[0]->E[iin+1]*dir)
        tritype = 1;
        else
        tritype = 2;
      /*} */

      /* make two triangles */
      if (tritype == 1) {
        t->tris[tri].p[0] = ii;
        t->tris[tri].p[1] = iin;
        t->tris[tri].p[2] = ii+1;
        t->tris[tri+1].p[0] = ii+1;
        t->tris[tri+1].p[1] = iin;
        t->tris[tri+1].p[2] = iin+1;
      }
      else {
        t->tris[tri].p[0] = ii;
        t->tris[tri].p[1] = iin;
        t->tris[tri].p[2] = iin+1;
        t->tris[tri+1].p[0] = ii;
        t->tris[tri+1].p[1] = iin+1;
        t->tris[tri+1].p[2] = ii+1;
      }

      /* ensure anti-clockwise order */
      if (dir < 0) {
        tmp = t->tris[tri].p[1];
        t->tris[tri].p[1] = t->tris[tri].p[2];
        t->tris[tri].p[2] = tmp;
        tmp = t->tris[tri+1].p[1];
        t->tris[tri+1].p[1] = t->tris[tri+1].p[2];
        t->tris[tri+1].p[2] = tmp;
      }

      /* go to next tris */
      tri += 2;
    }
  }

  return 0;
}
