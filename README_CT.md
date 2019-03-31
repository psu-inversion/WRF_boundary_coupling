  <DESCRIPTION>
  
#  Program wrfchembc.
  
For more information, contact Martha Butler or Thomas Lauvaux (tul5@psu.edu), The Pennsylvania State University.

Original source of unknown date. Authors: Rainer Schmitz (University of Chile - Santiago, Chile) Steven Peckham (NOAA/ESRL/GSD - Boulder, CO)
  
   This code has been liberally adapted at Penn State to support the 'nesting' of 
   a WRF domain within a global model, for the purpose of populating the boundaries
   (in wrfbdy) with a tracer from the global model.
   wrfchembc has versions specific to the global model being used for the problem.
   wrfchembc is executed once for each tracer to be populated.
   Code can also be used to initialize tracers with a constant value
         (options now are 300 ppm or 0 ppm, assuming the tracer is CO2)
         Use the 300 ppm option when the tracer will be used for surface fluxes
         with both source and sink values to avoid 'losing' tracer.
  
   Assumes existence of a wrfinput file and matching wrfbdy file populated with
     meteorological values from real.exe with both of these files having the
     tracers already defined. (It is possible to add the tracers to the
     wrfbdy file if they do not already exist, but this is not recommended.) 
  
   This is the main program of WRF/Chem global boundary condition generation code, and is
   specific to the global model being used.  The tracer for each execution is provided in 
   the 'control' namelist.  This code calls procedures in the module for the global model
   and in the separate module dealing with wrf procedures. 
   Only the wrfbdy file tracers are updated in this code. 
  
   The global file consists of a two-month concatenation of source data.
   Sourcing begins at the time-step within the 2-month file that corresponds
   to the beginning date of the wrfinput/wrfbdy files.
   It is up to the user to specify the correct 2-month global source file
   in the namelist control.
   For the TM5 and GEOS-Chem based files, the global file time increment (`dt_g`) is 3 hours.
   The wrfbdy time increment (`dt_w`) is 6 hours.
   Appropriate global file time incrments will be averaged to the wrfbdy time step.
   This code can accommodate either 00 UTC or 12 UTC start times.
  
   The options for global models at this time are:
          CarbonTracker   (other TM5 based global models, such as those from Sourish Basu,
                            will use this as a model, with some modifications)
          GEOS-Chem       (a CMS-Flux system version, which will also be used as a model for
                             files from Andrew Schuh, with some modifications)
          PCTM            (a PCTM system version, which requires re-orienting PCTM files
                             to surface-top before use) 
 
  The wrfchembc package consists of this 'main' code plus:
    `module_CT_lib.f90`:       but, use the module specific to global model
                             currently either `module_GEOS_lib.f90` or `module_CT_lib.f90`
                             Most procedures in the global model module will be referred to as 'global'
    `module_wrfchem_lib.f90`:  This contains wrf specific procedures and should not need to be modified.
    `wrfchembc_namelist.input`
    `Makefile`
  
   The namelist input data file for the program is structured as follows with examples:

```fortran
   &control
   
   dir_wrf = '/data/wrfchem_data/'                   ! Directory containing the wrfinput and wrfbdy files
   fnb_wrf  = 'wrfbdy_d01'                           ! WRF boundary condition data file
   fni_wrf  = 'wrfinput_d01'                         ! WRF initial condition data file
  
   chem_bc_opt = 1 (or 2 or 3)                       ! 1=global model, 2=fixed constant, 3=zero
  
   specname = 'tracer_8'                             
   
   !identification of the tracer to be updated in wrfbdy
   !specname = 'tracer_n' where n is the tracer number
   !tracer name must match the variable name in wrfbdy file
  
   The following are used only for chem_bc_opt = 1
   and may be omitted for options 2 and 3
   
   dir_global = '/data/global_data/'                 # Global model data directory example
   fn_global  = 'CT2016_2010-0607.nc'                # Global model data file name example
```
   
   To compile, type `make`. The make file is set-up to compile the code for a 
   single processor using the Intel compiler (We assume you are using a 
   linux cluster).  The Makefile will need to be modified to work with other 
   computer systems.
   Makefile in this directory is specific to PSU Meteorology Linux cluster use
  
   To execute the program, type 
  
   `wrfchembc_CT < wrfchembc_namelist.input`
  
 

