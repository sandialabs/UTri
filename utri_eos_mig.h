/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of UTri, substitue the full text  */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef UTRI_EOS_MIG_H
#define UTRI_EOS_MIG_H

#include "utri_FC.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
 * Note that all functions in the UTri EOS MIG interface assume
 * FORTRAN linking, and as such all variables are passed by reference.
 *
 */

/*
 * Initialize an instance of the utri eos library and return an
 * integer handle for it.
 *
 */
int utri_eos_handle_open();

/*
 * Free resources associated with a utri eos handle.
 *
 */
void utri_eos_handle_close( const int handle );

/*
 * get a string describing the passed error status
 *
 * inputs:
 *   error_code -- pointer to integer containing the error code. Code
 *                 should be one from the utri_error_codes enum.
 *
 * return:
 *   A null terminated character string containing the error
 *   description is returned. If non-null, the memory must be freed by
 *   the caller.
 *
 */
char * utri_eos_status( const int * error_code );

/*
 * The initialization function for a tabular model. Upon successful
 * return, the ui and dc arrays will define a tabular model that may
 * be passed to the other utri_eos calls. All arrays must be allocated
 * by the caller.
 *
 * inputs:
 *   ui -- double array of input parameters (refer to ui_index_names enum for size)
 *   gc -- double array (unused)
 *   dc -- double array (refer to dc_index_names for size)
 *         Upon entry the first seven entries contain factors for converting from
 *         the host code units into SI. For example, if the host is using CGS units
 *         then one would set dc[0] = 100.
 *         dc[0] => length unit
 *         dc[1] => mass unit
 *         dc[2] => time unit
 *         dc[3] => temperature unit
 *         dc[4] => discrete amount unit
 *         dc[5] => electric current unit
 *         dc[6] => luminosity unit
 *   uc -- constant character array containing the file name for the desired table.
 *         This file name will have the prefix '-xx.nc' appended to it where 'xx' is
 *         replaced by the abbreviation for the specified independent variable space.
 *   mdc -- integer giving the available size in the dc array
 *
 * outputs:
 *   dc  -- double array containing the model's derived constants (refer to dc_index_names)
 *   ndc -- integer giving the amount of space used in the dc array
 *   vi  -- double array of length 5 containing initial values for the model. Depending
 *          on the material, these may not all be sensible:
 *          vi[0] == initial density
 *          vi[1] == initial temperature
 *          vi[2] == minimum solid density
 *          vi[3] == melting temperature
 *          vi[4] == melting density
 */
void utri_eos_init( double * ui,
                    double * gc,
                    double * dc,
                    const char * uc,
                    int * mdc,
                    int * ndc,
                    double * vi,
                    int * handle,
                    int * error_code );

/*
 * Cleanup function for a utri tabular model. Upon return, the model
 * described by the ui and dc arrays will have been completely freed,
 * invalidating those arrays for use in calls to the other utri_eos
 * functions except for utri_eos_init.
 *
 * inputs:
 *   ui -- double array of model parameters
 *   gc -- double array (unused)
 *   dc -- double array of model derived constants
 *
 */  
void utri_eos_delete( double * ui,
                      double * gc,
                      double * dc );

/*
 * The pressure-{enthalpy,entropy,internal energy} space interpolation
 * routine. This routine performs vectorized look ups on an
 * unstructured triangular data table. Upon return the interpolated
 * quantities are found in the scratch array.
 *
 * inputs:
 *   mc -- integer vectorization stride length
 *   nc -- integer number of points to evaluate in the stride (requires nc <= mc)
 *   ui -- double array of input parameters (as from after a call to utri_eos_init)
 *   dc -- double array of derived constants (as from after call to utri_eos_init)
 *   pressure -- double array (length mc) of input pressures
 *   xvar     -- double array (length mc) of input second variable
 *   compute_vars -- bit wise or of desired variables to compute (see utri_dvar_names enum)
 *   int_scratch -- integer array (length 2*mc) used for interpolation calculations
 *
 * outputs:
 *   scratch     -- double array (length 19*mc) used for interpolation calculations.
 *                  Returns interpolated quantities with the following starting indices:
 *                  scratch[0*mc] => enthalpy
 *                  scratch[1*mc] => entropy
 *                  scratch[2*mc] => density
 *                  scratch[3*mc] => temperature
 *                  scratch[4*mc] => Gibbs free energy
 *                  scratch[5*mc] => Helmholtz free energy
 *                  scratch[6*mc] => adiabatic sound speed
 *                  scratch[7*mc] => isothermal bulk modulus
 *                  scratch[8*mc] => isobaric heat capacity
 *                  scratch[9*mc] => isothermal heat capacity
 *                  scratch[10*mc] => viscosity
 *                  scratch[11*mc] => thermal conductivity
 *                  scratch[12*mc] => phase
 *   error -- integer containing an error code. Non-zero indicates an
 *            error, see utri_error_codes enum.
 *
 */  
#define UTRI_PXDRV_SCRATCH_SIZE 19
#define UTRI_PXDRV_ISCRATCH_SIZE 2

void utri_eos_pxdrv( int * mc,
                     int * nc,
                     double * ui,
                     double * dc,
                     double * pressure,
                     double * xvar,
                     int * compute_vars,
                     double * scratch,
                     int * int_scratch,
                     int * error );

/*
 * Routine to obtain saturation values at a given pressure or
 * temperature along a given phase boundary. Upon return the values
 * for each phase are given in sat1 and sat2.
 *
 * input:
 *   mc -- integer vectorization stride length
 *   nc -- integer number of points to evaluate in the stride (requires nc <= mc)
 *   ui -- double array of input parameters (as from after a call to utri_eos_init)
 *   dc -- double array of derived constants (as from after call to utri_eos_init)
 *   ptvar    -- double array of pressures or temperatures at which to evaluate the boundary
 *   pttype   -- integer denoting whether ptvars are pressures (==0) or temperatures (==1)
 *   boundary -- integer denoting the boundary on which to evaluate
 *   compute_vars -- bit wise or of desired variables to compute (see utri_dvar_names enum)
 *   int_scratch -- integer array (length 2*mc) used for interpolation calculations
 *   scratch     -- double array (length 5*mc) used for interpolation calculations.
 *
 * outputs:
 *   sat1 -- double array (length 13*mc) containing interpolated saturation quantities
 *           for the first phase at the following starting indices:
 *           0*mc => enthalpy
 *           1*mc => entropy
 *           2*mc => density
 *           3*mc => temperature
 *           4*mc => Gibbs free energy
 *           5*mc => Helmholtz free energy
 *           6*mc => adiabatic sound speed
 *           7*mc => isothermal bulk modulus
 *           8*mc => isobaric heat capacity
 *           9*mc => isothermal heat capacity
 *           10*mc => viscosity
 *           11*mc => thermal conductivity
 *           12*mc => phase
 *   sat2 -- double array (length 13*mc) containing interpolated saturation quantities
 *           for the second phase with the same indices as sat1.
 *
 */  
#define UTRI_PTSAT_DATA_SIZE 13
#define UTRI_PTSAT_SCRATCH_SIZE 5
#define UTRI_PTSAT_ISCRATCH_SIZE 2

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
                     int * error );

/*
 * The density-temperature space interpolation routine. This routine
 * performs vectorized look ups on an unstructured triangular data
 * table.
 *
 * inputs:
 *   mc -- integer vectorization stride length
 *   nc -- integer number of points to evaluate in the stride (requires nc <= mc)
 *   ui -- double array of input parameters (as from after a call to utri_eos_init)
 *   gc -- double array (unused)
 *   dc -- double array of derived constants (as from after call to utri_eos_init)
 *   density     -- double array (length mc) of input densities
 *   temperature -- double array (length mc) of input temperatures
 *   alpha       -- double array (length mc*num_extra) of input extra porosity variables (unused)
 *   num_extra   -- integer number of extra variables in the alpha array
 *   scratch     -- double array (length 14*mc) used for interpolation calculations
 *   int_scratch -- integer array (length 2*mc) used for interpolation calculations
 *
 * outputs:
 *   pressure    -- double array (length mc) of output pressures
 *   energy      -- double array (length mc) of output internal energies
 *   sound_speed -- double array (length mc) of output adiabatic sound speeds
 *   scratch     -- Returns interpolated quantities with the following starting indices:
 *                  scratch[0*mc] => dP/dR
 *                  scratch[1*mc] => dP/dT
 *                  scratch[2*mc] => dE/dT
 *                  scratch[3*mc] => dE/dR
 *   error       -- integer pointer to error code
 */  
#define UTRI_RTDRV_SCRATCH_SIZE 14
#define UTRI_RTDRV_ISCRATCH_SIZE 2

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
                     int * error );

/*
 * The density-internal energy space interpolation routine. This
 * routine performs vectorized look ups on an unstructured triangular
 * data table.
 *
 * inputs:
 *   mc -- integer vectorization stride length
 *   nc -- integer number of points to evaluate in the stride (requires nc <= mc)
 *   ui -- double array of input parameters (as from after a call to utri_eos_init)
 *   gc -- double array (unused)
 *   dc -- double array of derived constants (as from after call to utri_eos_init)
 *   density     -- double array (length mc) of input densities
 *   energy      -- double array (length mc) of input internal energies
 *   alpha       -- double array (length mc*num_extra) of input extra porosity variables (unused)
 *   num_extra   -- integer number of extra variables in the alpha array
 *   pressure_flag -- integer flag to skip pressure calculation (==1) otherwise zero
 *   scratch     -- double array (length 14*mc) used for interpolation calculations
 *   int_scratch -- integer array (length 2*mc) used for interpolation calculations
 *
 * outputs:
 *   pressure    -- double array (length mc) of output pressures
 *   temperature -- double array (length mc) of output temperatures
 *   sound_speed -- double array (length mc) of output adiabatic sound speeds
 *   scratch     -- Returns interpolated quantities with the following starting indices:
 *                  scratch[0*mc] => dP/dR
 *                  scratch[1*mc] => dP/dT
 *                  scratch[2*mc] => dE/dT
 *                  scratch[3*mc] => dE/dR
 *   error       -- integer pointer to error code
 *
 */  
#define UTRI_REDRV_SCRATCH_SIZE 14
#define UTRI_REDRV_ISCRATCH_SIZE 2

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
                     int * error );

/*
 * This is a convenience enum to make the ui array indices easier to understand.
 *
 * The description of the parameters is given below
 */
enum ui_index_names {
  IVAR, /* requested independent variable space (see utri_ivar_names enum) */
  DVAR, /* requested dependent variables (see utri_dvar_names enum) */
  SR,   /* density scaling factor */
  R0,   /* initial density (can be set from table) */
  T0,   /* initial temperature (can be set from table) */
  RMIN, /* minimum solid density */
  Z,    /* average atomic number (set from table) */
  A,    /* average atomic weight (set from table) */
  ESFT, /* energy shift for model (can be calculated automatically */
  TYP,  /* model type (TYP==2 only supported value) */
  NXI,  /* number of desired table modes */
  XI0,  /* coefficient for table mode 0 */
  XI1,  /* coefficient for table mode 1 */
  XI2,  /* coefficient for table mode 2 */
  XI3,  /* coefficient for table mode 3 */
  XI4,  /* coefficient for table mode 4 */
  XI5,  /* coefficient for table mode 5 */
  XI6,  /* coefficient for table mode 6 */
  XI7,  /* coefficient for table mode 7 */
  XI8,  /* coefficient for table mode 8 */
  XI9,  /* coefficient for table mode 9 */
  NUM_UI_INDEX_NAMES
};

/*
 * This is a convenience enum to make the dc array indices easier to understand.
 *
 * Note that codes should not modify any values in the dc array as they are
 * set by the init routine. However, the NUM_DC_INDEX_NAMES value may be used
 * to allocate a properly sized array.
 *
 */
enum dc_index_names {
  KEY,
  CV0,
  CS0,
  RCNV,
  TCNV,
  PCNV,
  ECNV,
  CSCNV,
  MUCNV,
  LCNV,
  RCNVI,
  TCNVI,
  PCNVI,
  ECNVI,
  SQRT_SR,
  NUM_DC_INDEX_NAMES
};

/*
 * This enum gives the supported independent variable types
 *
 */
enum utri_ivar_names {
  I_RT,  /* density-temperature mesh */
  I_RE,  /* density-internal energy mesh */
  I_PH,  /* pressure-enthalpy mesh */
  I_PS,  /* pressure-entropy mesh */
  I_PE,  /* pressure-internal energy mesh */
  NUM_UTRI_IVAR_NAMES
};

/*
 * This enum gives the supported dependent variable types
 *
 */
enum utri_dvar_names {
  V_R  =   0x1,  /* density */
  V_T  =   0x2,  /* temperature */
  V_P  =   0x4,  /* pressure */
  V_G  =   0x8,  /* Gibb's free energy */
  V_F  =  0x10,  /* Helmholtz free energy */
  V_H  =  0x20,  /* enthalpy */
  V_E  =  0x40,  /* internal energy */
  V_S  =  0x80,  /* entropy */
  V_CP = 0x100,  /* heat capacity at constant pressure */
  V_CV = 0x200,  /* heat capacity at constant volume */
  V_CS = 0x400,  /* adiabatic sound speed */
  V_KT = 0x800   /* isothermal bulk modulus */
};

/*
 * This enum gives the location of dependent variable types in the
 * scratch array returned by the pxdrv routine and the sat[12] arrays
 * returned by the ptsat routine.
 *
 */
enum utri_px_dat_locs {
  PXD_H,  /* enthalpy                  */
  PXD_S,  /* entropy                   */
  PXD_R,  /* density                   */
  PXD_T,  /* temperature               */
  PXD_G,  /* Gibbs free energy         */
  PXD_F,  /* Helmholtz free energy     */
  PXD_CS, /* adiabatic sound speed     */
  PXD_KT, /* isothermal bulk modulus   */
  PXD_CP, /* isobaric heat capacity    */
  PXD_CV, /* isothermal heat capacity  */
  PXD_M,  /* viscosity                 */
  PXD_L,  /* thermal conductivity      */
  PXD_PH, /* phase                     */
  NUM_UTRI_PX_DAT_LOCS
};

/*
 * This enum gives error codes
 *
 */
enum utri_error_codes {
  /* no error */
  UTRI_SUCCESS           = 0x00,
  /* out of table bounds clip warnings */
  UTRI_ERROR_CLIP_X_LOW  = 0x01,
  UTRI_ERROR_CLIP_X_HIGH = 0x02,
  UTRI_ERROR_CLIP_Y_LOW  = 0x04,
  UTRI_ERROR_CLIP_Y_HIGH = 0x08,
  /* general error type */
  UTRI_FATAL             = 0x10,
  UTRI_WARNING           = 0x20,
  UTRI_INFO              = 0x40,
  /* end point */
  UTRI_ERROR_CODES_SIZE  = 0x40
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
