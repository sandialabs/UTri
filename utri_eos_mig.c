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
#include <math.h>

#include "utri_eos_support.h"
#include "utri_eos_mig.h"

/*
 * Initialize an instance of the utri eos library and return an
 * integer handle for it.
 *
 */
int utri_eos_handle_open()
{
  return utri_eos_handle_init(NULL,NULL,1);
}

/*
 * Close and free resources associated with a utri eos handle.
 *
 */
void utri_eos_handle_close( const int handle )
{
  utri_eos_handle_free(handle);
}

/*
 * This struct and its static instantiation holds the tables that have
 * been read into memory. They are identified by the keys which is the
 * name of the data file.  The model initialization looks up the
 * passed key and if not found reads the new table in from disk. The
 * key's integer is saved for fast look up in driver routines.  As
 * such, when allocating space for a new table, ordering of the keys
 * and tables must not be altered.
 */
struct utri_eos_table_data
{
  int n;
  eos_table_utri_t *tables;
};

static struct utri_eos_table_data utri_eos_tables = { 0, NULL };

/* user must free char* memory */
char * utri_eos_status( const int * error_code )
{
  /* concatenate strings corresponding to error code bit fields */
  int ec = *error_code;
  if (ec & UTRI_ERROR_CLIP_X_LOW)
    reportUTriError("UTRI Error: Clip X Low");
  if (ec & UTRI_ERROR_CLIP_X_HIGH)
    reportUTriError("UTRI Error: Clip X High");
  if (ec & UTRI_ERROR_CLIP_Y_LOW)
    reportUTriError("UTRI Error: Clip Y Low");
  if (ec & UTRI_ERROR_CLIP_Y_HIGH)
    reportUTriError("UTRI Error: Clip Y High");

  /* for general status, only report the most severe */
  if      (ec & UTRI_FATAL)   reportUTriError("UTRI Fatal:");
  else if (ec & UTRI_WARNING) reportUTriError("UTRI Warning:");
  else if (ec & UTRI_INFO)    reportUTriError("UTRI Info:");
  
  return utri_error_string_get();
}

/*
 * A helper function to get the file names out of the passed uc array
 *
 * The base file name is passed. The full file name is dependent upon
 * the independent variable type ivtype:
 * ivtype == I_R{T,E}: for rx tables, respectively append "-rx.nc"
 * ivtype == I_P{H,S,E}: for px tables, respectively append "-px.nc"
 *
 */
int utri_eos_uc_strings( const char * uc,
                         const int ivtype,
                         char ** file )
{
  int i,length;
  int fail;

  /*
   * as a hack to support a fortran string being passed, possibly
   * tokenized with '|', we assume that the string ends either when
   * finding a null, blank, or '|' character
   */
  length = 0;
  i = 0;
  while (1) {
    if (uc[i] == '\0' || uc[i] == ' ' || uc[i] == '|') {
      length = i;
      break;
    }
    i++;
  }

  /* check if no name found */
  if (length < 1) return -1;

  fail = 0;
  if ((*file = (char *)malloc(sizeof(char)*(length+7))) == NULL) fail = 1;
  /* memory allocation failure */
  if (fail) return -2;

  strncpy(*file,uc,length);
  if      (ivtype == I_RT) strcpy(&(*file)[length],"-rt.nc");
  else if (ivtype == I_RE) strcpy(&(*file)[length],"-re.nc");
  else if (ivtype == I_PH) strcpy(&(*file)[length],"-ph.nc");
  else if (ivtype == I_PS) strcpy(&(*file)[length],"-ps.nc");
  else if (ivtype == I_PE) strcpy(&(*file)[length],"-pe.nc");
  else return -3; /* invalid ivtype */

  return 0;
}

/*
 * The internal state variable initialization routine
 *
 */
void utri_eos_isvinit( double * ui,
                       double * gc,
                       double * dc,
                       int * num_extra,
                       char * namea,
                       char * keya,
                       double * rinit,
                       double * rdim,
                       int * iadvct,
                       int * itype )
{
  *num_extra = 0;
}

/*
 * The density-temperature based driver routine
 *
 * scratch is mc*14 in length
 * int_scratch is mc*2 in length
 *
 */
void utri_eos_rtdrv( int * mc,
                     int * nc,
                     double * ui,
                     double * gc,
                     double * dc,
                     double * density,
                     double * temperature,
                     double * alpha,
                     int * num_extra,
                     double * pressure,
                     double * energy,
                     double * sound_speed,
                     double * scratch,
                     int * int_scratch,
                     int * error )
{
  int table;
  int i;
  int n,stride;
  eos_data_points_t f;
  int clip[4];
  int logvars;

  n = *nc;
  stride = *mc;
  *error = UTRI_SUCCESS;

  /* table key */
  table = dc[KEY];

  /* independent vars in logscale? */
  logvars = utri_eos_tables.tables[table].flags[0];

  /* adjust density for units and scaling and temperature for units */
  for (i=0;i<n;i++) {
    scratch[i+7*stride] = density[i]*dc[RCNVI]*ui[SR];
    scratch[i+8*stride] = temperature[i]*dc[TCNVI];
  }

  /* convert to log space if necessary */
  if (logvars == 1) {
    for (i=0;i<n;i++) {
      scratch[i+7*stride] = log(scratch[i+7*stride]);
      scratch[i+8*stride] = log(scratch[i+8*stride]);
    }
  }

  /* fill eos_data_points structure */
  f.n = stride;
  f.phase = &int_scratch[0];
  f.E = energy;
  f.p = pressure;
  f.cs = sound_speed;
  f.dpdr = &scratch[0*stride];
  f.dpdT = &scratch[1*stride];
  f.dEdT = &scratch[2*stride];
  f.dEdr = &scratch[3*stride];
  f.dcsdr = &scratch[4*stride];
  f.visc  = &scratch[5*stride];
  f.thcon = &scratch[6*stride];

  /* zero clip counters */
  for (i=0;i<4;i++) clip[i] = 0;

  eostableutrievaluatec0(&utri_eos_tables.tables[table],stride,n,&scratch[7*stride],&scratch[8*stride],&f,1,clip,&scratch[9*stride],&int_scratch[stride]);

  /* issue clipping warnings */
  if (clip[0] > 0) *error |= UTRI_ERROR_CLIP_X_LOW;
  if (clip[1] > 0) *error |= UTRI_ERROR_CLIP_X_HIGH;
  if (clip[2] > 0) *error |= UTRI_ERROR_CLIP_Y_LOW;
  if (clip[3] > 0) *error |= UTRI_ERROR_CLIP_Y_HIGH;

  /* scale output quantities and convert back to code units */
  for (i=0;i<n;i++) {
    f.p[i]     *= dc[PCNV];
    f.E[i]      = f.E[i]*dc[ECNV]*ui[SR]+ui[ESFT];
    f.cs[i]    *= dc[CSCNV]*dc[SQRT_SR];
    f.dpdr[i]  *= dc[PCNV]*dc[RCNVI]*ui[SR];
    f.dpdT[i]  *= dc[PCNV]*dc[TCNVI];
    f.dEdT[i]  *= dc[ECNV]*dc[TCNVI]*ui[SR];
    f.dEdr[i]  *= dc[ECNV]*dc[RCNVI]*ui[SR]*ui[SR];
    f.dcsdr[i] *= dc[CSCNV]*dc[RCNVI]*dc[SQRT_SR]*ui[SR];
    f.visc[i]  *= dc[MUCNV];
    f.thcon[i] *= dc[LCNV]*ui[SR];
    scratch[i+7*stride] = f.phase[i];
  }
}
  
/*
 * The density-internal energy based driver routine
 *
 * scratch is mc*14 in length
 * int_scratch is mc*2 in length
 *
 */
void utri_eos_redrv( int * mc,
                     int * nc,
                     double * ui,
                     double * gc,
                     double * dc,
                     double * density,
                     double * energy,
                     double * alpha,
                     int * num_extra,
                     int * pressure_flag,
                     double * pressure,
                     double * temperature,
                     double * sound_speed,
                     double * scratch,
                     int * int_scratch,
                     int * error )
{
  int table;
  int i;
  int n,stride;
  eos_data_points_t f;
  int clip[4];
  int logvars;

  n = *nc;
  stride = *mc;
  *error = UTRI_SUCCESS;

  /* table key */
  table = dc[KEY];

  /* independent vars in logscale? */
  logvars = utri_eos_tables.tables[table].flags[0];

  /* adjust density for units and scaling and energy for units */
  for (i=0;i<n;i++) {
    scratch[i+7*stride] = density[i]*dc[RCNVI]*ui[SR];
    scratch[i+8*stride] = (energy[i]-ui[ESFT])*dc[ECNVI]/ui[SR];
  }

  /* convert to log space if necessary */
  if (logvars == 1) {
    for (i=0;i<n;i++) {
      scratch[i+7*stride] = log(scratch[i+7*stride]);
      scratch[i+8*stride] = log(scratch[i+8*stride]);
    }
  }

  /* fill eos_data_points structure */
  f.n = stride;
  f.phase = &int_scratch[0];
  f.T = temperature;
  f.p = pressure;
  f.cs = sound_speed;
  f.dpdr = &scratch[0*stride];
  f.dpdT = &scratch[1*stride];
  f.dEdT = &scratch[2*stride];
  f.dEdr = &scratch[3*stride];
  f.dcsdr = &scratch[4*stride];
  f.visc  = &scratch[5*stride];
  f.thcon = &scratch[6*stride];

  /* zero clip counters */
  for (i=0;i<4;i++) clip[i] = 0;

  eostableutrievaluatec0(&utri_eos_tables.tables[table],stride,n,&scratch[7*stride],&scratch[8*stride],&f,2,clip,&scratch[9*stride],&int_scratch[stride]);

  /* issue clipping warnings */
  if (clip[0] > 0) *error |= UTRI_ERROR_CLIP_X_LOW;
  if (clip[1] > 0) *error |= UTRI_ERROR_CLIP_X_HIGH;
  if (clip[2] > 0) *error |= UTRI_ERROR_CLIP_Y_LOW;
  if (clip[3] > 0) *error |= UTRI_ERROR_CLIP_Y_HIGH;
 
  /* scale output quantities and convert back to code units */
  for (i=0;i<n;i++) {
    f.p[i]     *= dc[PCNV];
    f.T[i]     *= dc[TCNV];
    f.cs[i]    *= dc[CSCNV]*dc[SQRT_SR];
    f.dpdr[i]  *= dc[PCNV]*dc[RCNVI]*ui[SR];
    f.dpdT[i]  *= dc[PCNV]*dc[TCNVI];
    f.dEdT[i]  *= dc[ECNV]*dc[TCNVI]*ui[SR];
    f.dEdr[i]  *= dc[ECNV]*dc[RCNVI]*ui[SR]*ui[SR];
    f.dcsdr[i] *= dc[CSCNV]*dc[RCNVI]*dc[SQRT_SR]*ui[SR];
    f.visc[i]  *= dc[MUCNV];
    f.thcon[i] *= dc[LCNV]*ui[SR];
    scratch[i+7*stride] = f.phase[i];
  }
}

/*
 * The pressure-{enthalpy/entropy/internal energy} based driver routine
 *
 * scratch is mc*19 in length
 * int_scratch is mc*2 in length
 *
 */
void utri_eos_pxdrv( int * mc,
                     int * nc,
                     double * ui,
                     double * dc,
                     double * pressure,
                     double * xvar,
                     int * compute_vars,
                     double * scratch,
                     int * int_scratch,
                     int * error )
{
  int table;
  int i;
  int n,stride;
  eos_data_points_t f;
  int clip[4];
  int ivtype;
  int logvars;

  n = *nc;
  stride = *mc;
  ivtype = ui[IVAR];
  *error = UTRI_SUCCESS;

  /* table key */
  table = dc[KEY];

  /* independent vars in logscale? */
  logvars = utri_eos_tables.tables[table].flags[0];

  /* adjust pressure for units and xvar for scaling and units */
  for (i=0;i<n;i++) {
    scratch[i+12*stride] = pressure[i]*dc[PCNVI];
    if (ivtype == I_PH)      scratch[i+13*stride] = (xvar[i]-ui[ESFT])*dc[ECNVI]/ui[SR];
    else if (ivtype == I_PS) scratch[i+13*stride] = xvar[i]*dc[ECNVI]*dc[TCNV]/ui[SR];
    else                     scratch[i+13*stride] = (xvar[i]-ui[ESFT])*dc[ECNVI]/ui[SR];
  }

  /* convert to log space if necessary */
  if (logvars == 1) {
    for (i=0;i<n;i++) {
      scratch[i+12*stride] = log(scratch[i+12*stride]);
      scratch[i+13*stride] = log(scratch[i+13*stride]);
    }
  }

  /* fill eos_data_points structure */
  f.n = stride;
  f.phase = &int_scratch[0];
  if (ivtype == I_PH) {
    f.E = &scratch[0*stride];
    f.S = &scratch[1*stride];
  }
  else if (ivtype == I_PS) {
    f.E = &scratch[0*stride];
    f.dcsdr = &scratch[1*stride]; /* really H*/
  }
  else {
    f.dcsdr = &scratch[0*stride]; /* really H*/
    f.S = &scratch[1*stride];
  }
  f.rho = &scratch[2*stride];
  f.T = &scratch[3*stride];
  f.dEdr = &scratch[4*stride]; /* really G*/
  f.F = &scratch[5*stride];
  f.cs = &scratch[6*stride];
  f.dpdr = &scratch[7*stride];
  f.dpdT = &scratch[8*stride];
  f.dEdT = &scratch[9*stride];
  f.visc = &scratch[10*stride];
  f.thcon = &scratch[11*stride];

  /* zero clip counters */
  for (i=0;i<4;i++) clip[i] = 0;

  eostableutrievaluatec0(&utri_eos_tables.tables[table],stride,n,&scratch[12*stride],&scratch[13*stride],&f,ivtype+2,clip,&scratch[14*stride],&int_scratch[stride]);

  /* issue clipping warnings */
  if (clip[0] > 0) *error |= UTRI_ERROR_CLIP_X_LOW;
  if (clip[1] > 0) *error |= UTRI_ERROR_CLIP_X_HIGH;
  if (clip[2] > 0) *error |= UTRI_ERROR_CLIP_Y_LOW;
  if (clip[3] > 0) *error |= UTRI_ERROR_CLIP_Y_HIGH;
 
  /* scale output quantities and convert back to code units */
  for (i=0;i<n;i++) {
    if (ivtype == I_PH) {
      f.E[i]  = f.E[i]*dc[ECNV]*ui[SR]+ui[ESFT];
      f.S[i] *= dc[ECNV]*dc[TCNVI];
    }
    else if (ivtype == I_PS) {
      f.E[i]  = f.E[i]*dc[ECNV]*ui[SR]+ui[ESFT];
      f.dcsdr[i]  = f.dcsdr[i]*dc[ECNV]*ui[SR]+ui[ESFT]; /* really H*/
    }
    else {
      f.dcsdr[i]  = f.dcsdr[i]*dc[ECNV]*ui[SR]+ui[ESFT]; /* really H*/
      f.S[i] *= dc[ECNV]*dc[TCNVI];
    }
    f.rho[i]  *= dc[RCNV]/ui[SR];
    f.T[i]    *= dc[TCNV];
    f.dEdr[i]     = f.dEdr[i]*dc[ECNV]*ui[SR]+ui[ESFT]; /* really G*/
    f.F[i]     = f.F[i]*dc[ECNV]*ui[SR]+ui[ESFT];
    f.cs[i]   *= dc[CSCNV]*dc[SQRT_SR];
    f.dpdr[i] *= dc[PCNV]; /* actually K_T */
    f.dpdT[i] *= dc[ECNV]*dc[TCNVI]*ui[SR]; /* actually C_P */
    f.dEdT[i] *= dc[ECNV]*dc[TCNVI]*ui[SR];
    f.visc[i]  *= dc[MUCNV];
    f.thcon[i] *= dc[LCNV]*ui[SR];
    scratch[i+12*stride] = f.phase[i];
  }
}

/*
 * The internal state variable quadrature routine
 *
 */
void utri_eos_isvquad( int * mc,
                       int * nc,
                       double * ui,
                       double * gc,
                       double * dc,
                       double * dt,
                       int * num_extra,
                       double * alph,
                       double * rho_0,
                       double * energy_0,
                       double * pressure_0,
                       double * temperature_0,
                       double * volume_change,
                       double * energy_change,
                       double * scratch )
{

}

/*
 * The cleanup routine which removes allocated tables from the table
 * array
 *
 */
void utri_eos_delete( double * ui,
                      double * gc,
                      double * dc )
{
  int table;
  int i;

  table = dc[KEY];
  eos_table_utri_delete(&utri_eos_tables.tables[table]);
  utri_eos_tables.tables[table].defaults = NULL;

  /* delete array if all entries empty */
  table = 0;
  for (i=0;i<utri_eos_tables.n;i++)
    if (utri_eos_tables.tables[i].defaults != NULL) table++;

  if (table == 0) {
    if (utri_eos_tables.tables) free(utri_eos_tables.tables);
    utri_eos_tables.tables = NULL;
    utri_eos_tables.n = 0;
  }

}

/*
 * The initialization routine which parses input and reads in the
 * needed data files
 *
 */
void utri_eos_init( double * ui,
                    double * gc,
                    double * dc,
                    const char * uc,
                    int * mdc,
                    int * ndc,
                    double * vi,
                    int * handle,
                    int * error_code )
{
  int status;
  int i,j;
  char *eosfile=NULL;
  int eoskey;
  int new;
  double rho_convert,t_convert,p_convert,e_convert,cs_convert,mu_convert,l_convert;
  int one=1,zero=0;
  double xv;
  double pi,ei,csi;
  double scratch[14];
  int int_scratch[2];
  char mesg[160];
  eos_table_utri_t * eostable = NULL;
  int ivtype;
  int logvars;
  eos_data_points_t * m0;
  int err;

  /* clear error status */
  *error_code = UTRI_SUCCESS;

  /*
   * Check for enough room in dc array
   *
   */
  *ndc = NUM_DC_INDEX_NAMES;
  if (*mdc < NUM_DC_INDEX_NAMES) {
    snprintf(mesg,160,"utri_eos_init: mdc = %d; needed %d",*mdc,*ndc);
    reportUTriError(mesg);
    *error_code |= UTRI_FATAL;
    return;
  }

  /*
   * Check for valid independent variable request
   *
   */
  ivtype = ui[IVAR];
  if (ivtype < 0 || ivtype >= NUM_UTRI_IVAR_NAMES) {
    snprintf(mesg,160,"utri_eos_init: invalid independent var type %d",ivtype);
    reportUTriError(mesg);
    *error_code |= UTRI_FATAL;
    return;
  }

  /*
   * Compute conversion factors
   *
   */
  if (ivtype == I_RT || ivtype == I_RE) {
    /*
     * RX space tables currently required to store info in SESAME units
     * (P in GPa, E in MJ/kg, rho in g/cm^3, T in K)
     */
    rho_convert = 1.e3*dc[1]/dc[0]/dc[0]/dc[0];
    t_convert = dc[3];
    p_convert = 1.e9*dc[1]/dc[0]/dc[2]/dc[2];
    e_convert = 1.e6*dc[0]*dc[0]/dc[2]/dc[2];
    cs_convert = 1.e3*dc[0]/dc[2];
    /* these two not used currently in this table type */
    mu_convert = 1.;
    l_convert = 1.;
  }
  else {
    /*
     * PX space tables currently required to store info in SI units
     *
     */
    rho_convert = dc[1]/dc[0]/dc[0]/dc[0];
    t_convert = dc[3];
    p_convert = dc[1]/dc[0]/dc[2]/dc[2];
    e_convert = dc[0]*dc[0]/dc[2]/dc[2];
    cs_convert = dc[0]/dc[2];
    mu_convert = dc[1]/dc[0]/dc[2];
    l_convert = dc[0]*dc[1]/dc[2]/dc[2]/dc[2]/dc[3];
  }

  dc[RCNV] = rho_convert;
  dc[TCNV] = t_convert;
  dc[PCNV] = p_convert;
  dc[ECNV] = e_convert;
  dc[CSCNV] = cs_convert;
  dc[MUCNV] = mu_convert;
  dc[LCNV] = l_convert;
  dc[RCNVI] = 1./rho_convert;
  dc[TCNVI] = 1./t_convert;
  dc[PCNVI] = 1./p_convert;
  dc[ECNVI] = 1./e_convert;

  /* check for sane density scaling */
  if (ui[SR] <= 0.) {
    snprintf(mesg,160,"utri_eos_init: non-positive density scaling factor %g reset to 1.0",ui[SR]);
    reportUTriError(mesg);
    *error_code |= UTRI_INFO;
    ui[SR] = 1.;
  }

  dc[SQRT_SR] = sqrt(ui[SR]);

  /* get file names from uc input */
  status = utri_eos_uc_strings(uc,ivtype,&eosfile);
  if (status != 0) {
    snprintf(mesg,160,"utri_eos_init: failed to generate input file names from base: %s\n",uc);
    reportUTriError(mesg);
    *error_code |= UTRI_FATAL;
    return;
  }

  /* find free key value */
  eoskey = -1;
  for (i=0;i<utri_eos_tables.n;i++)
    if (utri_eos_tables.tables[i].defaults == NULL)
      eoskey = i;
  
  /* reallocate data holders if needed */
  if (eoskey < 0) {
    eostable = (eos_table_utri_t *)realloc(utri_eos_tables.tables,(utri_eos_tables.n+1)*sizeof(eos_table_utri_t));
    if (eostable == NULL) {
      reportUTriError("utri_eos_init: failed to allocate memory for table");
      *error_code |= UTRI_FATAL;
      return;
    }
    utri_eos_tables.tables = eostable;
    eoskey = utri_eos_tables.n;
    utri_eos_tables.n++;
  }

  /* load table */
  status = eostableutrireadnc(eosfile,&utri_eos_tables.tables[eoskey],*handle);
  if (status != 0) {
    snprintf(mesg,160,"utri_eos_init: failed to read table file: %s\n",eosfile);
    utri_eos_tables.tables[eoskey].defaults = NULL;
    reportUTriError(mesg);
    *error_code |= UTRI_FATAL;
    return;
  }

  /* Check that the number of uncertain modes in the input table is at
     least as big as NXI */
  eostable = &utri_eos_tables.tables[eoskey];
  if (ui[NXI]+1 > eostable->nmodes ) {
    snprintf(mesg,160,"utri_eos_init: User requested too many uncertain modes which do not exist in the table %s\n",eosfile);
    reportUTriError(mesg);
    *error_code |= UTRI_FATAL;
    return;
  }

  /* mode 0 shortcut */
  m0 = eostable->modes[0];

  /* convert independent variables to log values if needed
   * this must be done before mode summation and rptree generation
   */
  logvars = eostable->flags[0];
  if (logvars) {
    if (ivtype == I_RT || ivtype == I_RE) {
      /* common variable */
      for (i=0;i<eostable->np;i++) m0->rho[i] = log(m0->rho[i]);
      /* unique variable */
      if (ivtype == I_RT) 
        for (i=0;i<eostable->np;i++) m0->T[i] = log(m0->T[i]);
      else if (ivtype == I_RE)
        for (i=0;i<eostable->np;i++) m0->E[i] = log(m0->E[i]);
    }
    else {
      /* common variable */
      for (i=0;i<eostable->np;i++) m0->p[i] = log(m0->p[i]);
      /* unique variable */
      if (ivtype == I_PH)
        for (i=0;i<eostable->np;i++) m0->H[i] = log(m0->H[i]);
      else if (ivtype == I_PS)
        for (i=0;i<eostable->np;i++) m0->S[i] = log(m0->S[i]);
      else if (ivtype == I_PE)
        for (i=0;i<eostable->np;i++) m0->E[i] = log(m0->E[i]);
    }
  }

  /* Modify the tables to account for uncertain coefficient found in
   * XI0,...,XI9.
   */
  for (i=0;i<ui[NXI];i++) {
    for (j=0;j<eostable->np;j++) {
      /* skip phase, it is invariant across the modes */
      if (logvars && (ivtype == I_RT || ivtype == I_RE))
        m0->rho[j]   += ui[NXI+i+1]*log(eostable->modes[i+1]->rho[j]);
      else 
        m0->rho[j]   += ui[NXI+i+1]*eostable->modes[i+1]->rho[j];
      if (logvars && ivtype == I_RT)
        m0->T[j]     += ui[NXI+i+1]*log(eostable->modes[i+1]->T[j]);
      else
        m0->T[j]     += ui[NXI+i+1]*eostable->modes[i+1]->T[j];
      if (logvars && (ivtype == I_PH || ivtype == I_PS || ivtype == I_PE))
        m0->p[j]     += ui[NXI+i+1]*log(eostable->modes[i+1]->p[j]);
      else
        m0->p[j]     += ui[NXI+i+1]*eostable->modes[i+1]->p[j];
      if (logvars && (ivtype == I_RE || ivtype == I_PE))
        m0->E[j]     += ui[NXI+i+1]*log(eostable->modes[i+1]->E[j]);
      else
        m0->E[j]     += ui[NXI+i+1]*eostable->modes[i+1]->E[j];
      if (logvars && ivtype == I_PS)
        m0->S[j]     += ui[NXI+i+1]*log(eostable->modes[i+1]->S[j]);
      else
        m0->S[j]     += ui[NXI+i+1]*eostable->modes[i+1]->S[j];
      m0->F[j]     += ui[NXI+i+1]*eostable->modes[i+1]->F[j];
      m0->dpdr[j]  += ui[NXI+i+1]*eostable->modes[i+1]->dpdr[j];
      m0->dpdT[j]  += ui[NXI+i+1]*eostable->modes[i+1]->dpdT[j];
      m0->dEdr[j]  += ui[NXI+i+1]*eostable->modes[i+1]->dEdr[j];
      m0->dEdT[j]  += ui[NXI+i+1]*eostable->modes[i+1]->dEdT[j];
      m0->cs[j]    += ui[NXI+i+1]*eostable->modes[i+1]->cs[j];
      if (logvars && ivtype == I_PH)
        m0->dcsdr[j] += ui[NXI+i+1]*log(eostable->modes[i+1]->dcsdr[j]);
      else
        m0->dcsdr[j] += ui[NXI+i+1]*eostable->modes[i+1]->dcsdr[j];
    }
  }

  if (ivtype == I_RT) {
    /* regenerate rptree for density-temperature table */
    if (calculate_rptree(eostable,0) != 0) {
      snprintf(mesg,160,"utri_eos_init: Failed to regenerate rptree for table: %s\n",eosfile);
      reportUTriError(mesg);
      *error_code |= UTRI_FATAL;
      return;
    }
  }
  else if (ivtype == I_RE) {
    /* regenerate rptree for density-energy table */
    if (construct_boundary_tris(eostable,2) != 0) {
      snprintf(mesg,160,"utri_eos_init: Failed to construct boundaries for table: %s\n",eosfile);
      reportUTriError(mesg);
      *error_code |= UTRI_FATAL;
      return;
    }

    if (calculate_rptree(eostable,1) != 0) {
      snprintf(mesg,160,"utri_eos_init: Failed to regenerate rptree for table: %s\n",eosfile);
      reportUTriError(mesg);
      *error_code |= UTRI_FATAL;
      return;
    }
  }
  else {
    /* regenerate rptree for pressure table */
    if (construct_boundary_tris(eostable,ivtype+2) != 0) {
      snprintf(mesg,160,"utri_eos_init: Failed to construct boundaries for table: %s\n",eosfile);
      reportUTriError(mesg);
      *error_code |= UTRI_FATAL;
      return;
    }

    if (calculate_rptree(eostable,ivtype+1) != 0) {
      snprintf(mesg,160,"eos_utri_init: Failed to regenerate rptree for table: %s\n",eosfile);
      reportUTriError(mesg);
      *error_code |= UTRI_FATAL;
      return;
    }
  }

  /* report clip boundaries */
  snprintf(mesg,80,"utri_eos_init: X boundaries %g to %g\n",
           (logvars==0)?utri_eos_tables.tables[eoskey].rptree.xbounds[0]:
           exp(utri_eos_tables.tables[eoskey].rptree.xbounds[0]),
           (logvars==0)?utri_eos_tables.tables[eoskey].rptree.xbounds[1]:
           exp(utri_eos_tables.tables[eoskey].rptree.xbounds[1]));
  reportUTriError(mesg);
  snprintf(mesg,80,"utri_eos_init: Y boundaries %g to %g\n",
           (logvars==0)?utri_eos_tables.tables[eoskey].rptree.ybounds[0]:
           exp(utri_eos_tables.tables[eoskey].rptree.ybounds[0]),
           (logvars==0)?utri_eos_tables.tables[eoskey].rptree.ybounds[1]:
           exp(utri_eos_tables.tables[eoskey].rptree.ybounds[1]));
  reportUTriError(mesg);
  *error_code |= UTRI_INFO;

  /* save key in dc array */
  dc[KEY] = eoskey;

  /* set/compute defaults and derived parameters */
  if (ui[R0] <= 0.) ui[R0] = utri_eos_tables.tables[eoskey].defaults[2]*dc[RCNV];
  if (ui[T0] <= 0.) ui[T0] = utri_eos_tables.tables[eoskey].defaults[3]*dc[TCNV];
  if (ui[Z] <= 0.) ui[Z] = utri_eos_tables.tables[eoskey].defaults[0];
  if (ui[A] <= 0.) ui[A] = utri_eos_tables.tables[eoskey].defaults[1];

  if (ivtype == I_RT) {
    xv = 0.;
    utri_eos_rtdrv(&one,&one,ui,gc,dc,&ui[R0],&ui[T0],&xv,&zero,&pi,&ei,&csi,scratch,int_scratch,&err);
    *error_code |= err;
    dc[CV0] = scratch[2];
    dc[CS0] = csi;

    /* default energy shift */
    new = rint(ui[TYP]);
    if (new != 2 && ui[ESFT] <= 0.) {
      ui[ESFT] = dc[CV0]*ui[T0];
      if (ei < 0.) ui[ESFT] -= ei;
    }
  }
  else {
    dc[CV0] = 0.;
    dc[CS0] = 0.;
  }

  /* set initial state defaults */
  vi[0] = ui[R0];
  vi[1] = ui[T0];
  vi[2] = 0.8*vi[0];/*ui[RMIN];*/
  vi[3] = 5.*vi[1];
  vi[4] = 0.8*vi[0];

  /* printf("vi %g %g %g %g %g\n",vi[0],vi[1],vi[2],vi[3],vi[4]); */

  /* free table strings */
  if (eosfile) free(eosfile);

}

/*
 * Routine to obtain saturated state values along a first order phase
 * transition line
 *
 * sat1 and sat2 are mc*11 in length
 * scratch is mc*5 in length
 * int_scratch is mc*2 in length
 *
 */
void utri_eos_ptsat( int * mc,
                     int * nc,
                     double * ui,
                     double * dc,
                     double * ptvar,
                     int * pttype,
                     int * boundary,
                     int * compute_vars,
                     double * sat1,
                     double * sat2,
                     double * scratch,
                     int * int_scratch,
                     int * error )
{
  reportUTriError("utri_eos_ptsat not yet implemented");
  *error = UTRI_FATAL;
}
