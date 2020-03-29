program main_bc_wrfchem_CT

  !<DESCRIPTION>
  !
  ! Program wrfchembc.
  !
  !     Original source of unknown date:
  ! Authors:  Rainer Schmitz (University of Chile - Santiago, Chile)
  !           Steven Peckham (NOAA/ESRL/GSD - Boulder, CO)
  !
  ! This code has been liberally adapted at Penn State to support the 'nesting' of 
  ! a WRF domain within a global model, for the purpose of populating the boundaries
  ! (in wrfbdy) with a tracer from the global model.
  ! wrfchembc has versions specific to the global model being used for the problem.
  ! wrfchembc is executed once for each tracer to be populated.
  ! Code can also be used to initialize tracers with a constant value
  !       (options now are 300 ppm or 0 ppm, assuming the tracer is CO2)
  !       Use the 300 ppm option when the tracer will be used for surface fluxes
  !       with both source and sink values to avoid 'losing' tracer.
  !
  ! Assumes existence of a wrfinput file and matching wrfbdy file populated with
  !   meteorological values from real.exe with both of these files having the
  !   tracers already defined. (It is possible to add the tracers to the
  !   wrfbdy file if they do not already exist, but this is not recommended.) 
  !
  ! This is the main program of WRF/Chem global boundary condition generation code, and is
  ! specific to the global model being used.  The tracer for each execution is provided in 
  ! the 'control' namelist.  This code calls procedures in the module for the global model
  ! and in the separate module dealing with wrf procedures. 
  ! Only the wrfbdy file tracers are updated in this code. 
  !
  ! The global file consists of a two-month concatenation of source data.
  ! Sourcing begins at the time-step within the 2-month file that corresponds
  ! to the beginning date of the wrfinput/wrfbdy files.
  ! It is up to the user to specify the correct 2-month global source file
  ! in the namelist control.
  ! For the TM5 and GEOS-Chem based files, the global file time increment (dt_g) is 3 hours.
  ! The wrfbdy time increment (dt_w) is 6 hours.
  ! Appropriate global file time incrments will be averaged to the wrfbdy time step.
  ! This code can accommodate either 00 UTC or 12 UTC start times.
  !
  !  The options for global models at this time are:
  !        CarbonTracker   (other TM5 based global models, such as those from Sourish Basu,
  !                          will use this as a model, with some modifications)
  !        GEOS-Chem       (a CMS-Flux system version, which will also be used as a model for
  !                           files from Andrew Schuh, with some modifications)
  !        PCTM            (a PCTM system version, which requires re-orienting PCTM files
  !                           to surface-top before use) 
 

  ! The wrfchembc package consists of this 'main' code plus:
  !  module_CT_lib.f90:       but, use the module specific to global model
  !                           currently either module_GEOS_lib.f90 or module_CT_lib.f90 
  !                           Most procedures in the global model module will be referred to as 'global'
  !  module_wrfchem_lib.f90:  This contains wrf specific procedures and should not need to be modified.
  !  wrfchembc_namelist.input
  !  Makefile
  !
  ! The namelist input data file for the program is structured as follows with examples:
  ! &control
  ! 
  ! dir_wrf = '/data/wrfchem_data/'                   ! Directory containing the wrfinput and wrfbdy files
  ! fnb_wrf  = 'wrfbdy_d01'                           ! WRF boundary condition data file
  ! fni_wrf  = 'wrfinput_d01'                         ! WRF initial condition data file
  !
  ! chem_bc_opt = 1 (or 2 or 3)                       ! 1=global model, 2=fixed constant, 3=zero
  !
  ! specname = 'tracer_8'                             ! identification of the tracer to be updated in wrfbdy
  !                                                   ! specname = 'tracer_n' where n is the tracer number
  !                                                   ! tracer name must match the variable name in wrfbdy file
  !
  !                                                   ! The following are used only for chem_bc_opt = 1
  !                                                   ! and may be omitted for options 2 and 3
  ! dir_global = '/data/global_data/'                 ! Global model data directory example
  ! fn_global  = 'CT2016_2010-0607.nc'                ! Global model data file name example
  ! 
  ! 
  ! /
  !
  ! To compile, type "make". The make file is set-up to compile the code for a 
  ! single processor using the Intel compiler (We assume you are using a 
  ! linux cluster).  The Makefile will need to be modified to work with other 
  ! computer systems.
  ! Makefile in this directory is specific to PSU Meteorology Linux cluster use
  !
  ! To execute the program, type 
  !
  !  wrfchembc_CT < wrfchembc_namelist.input
  !
  !</DESCRIPTION>

  use netcdf
  use module_wrfchem_lib
  use module_CT_lib        

  implicit none

  !     parameters - modify to suit
  integer, parameter :: maxsize = 100          ! for directories and filenames (and total length of dir + filename)
  integer, parameter :: litsize = 20           ! for tracer names
  integer, parameter :: dt_w = 6               ! the wrfbdy time increment in hours
  integer, parameter :: dt_g = 3               ! the global model source time increment in hours
  integer, parameter :: gincr = 8              ! number of global model time increments per day

  !     control variables (provided in namelist)
  character(len=maxsize) :: dir_global, dir_wrf
  character(len=maxsize) :: fnb_wrf, fni_wrf, fn_global
  integer                :: chem_bc_opt
  character(len=litsize) :: specname 

  namelist /control/ dir_wrf, fnb_wrf, fni_wrf, &
       chem_bc_opt, specname, dir_global, fn_global


  !     boundary conditions 
  real, allocatable :: bcxs(:,:,:)
  real, allocatable :: bcxe(:,:,:)
  real, allocatable :: bcys(:,:,:)
  real, allocatable :: bcye(:,:,:)
  !     boundary condition tendencies
  real, allocatable :: btxs(:,:,:)
  real, allocatable :: btxe(:,:,:)
  real, allocatable :: btys(:,:,:)
  real, allocatable :: btye(:,:,:)
  !     saved boundary conditions (for computing tendencies)
  real, allocatable :: bcxso(:,:,:)
  real, allocatable :: bcxeo(:,:,:)
  real, allocatable :: bcyso(:,:,:)
  real, allocatable :: bcyeo(:,:,:)


  !     other public working variables
  !     it_bdy is the 'time step' in the wrfbdy output file
  !     it_glob  is the 'time step' in the source global file (from file beginning)

  integer                :: it_bdy, it_glob
  character(len=maxsize) :: fnb, fni, fng    ! full names of files

  print *,"******************************************************"
  print *,"*   PROGRAM: WRF/CHEM GLOBAL BOUNDARIES              *"
  print *,"*              FOR WRF/CHEM VERSION 3.6.1            *"
  print *,"*                                                    *"
  print *,"*    PLEASE REPORT ANY BUGS TO us at PSU             *"
  print *,"*                                                    *"
  print *,"*             tul5@psu.edu                           *"
  print *,"*                                                    *"
  print *,"******************************************************"

  !     read control variables
  read( 5, nml=control )

  !     intialize module_wrfchem_lib
  fnb = trim(dir_wrf)//adjustl(fnb_wrf)      !wrfbdy
  fni = trim(dir_wrf)//adjustl(fni_wrf)      !wrfinput

  write(*,*) 'call to init_wrfchem_lib'
  !  open the wrf files, retrieving dimension data
  !      (will create boundary variables for this tracer in wrfbdy for this tracer if necessary)
  !   ntime (returned) is the total number of time increments in the wrfbdy file
  !   Other returned variables will be used in the interpolation routines
  call init_wrfchem_lib( fnb, fni, specname, ntime, ptop_wrf, znu, xlon, xlat )

  !     initialize global data module (if populating this tracer with global values)
  !    and compute horizontal and vertical interpolation coefficients 
  
  if(chem_bc_opt == 1) then
     fng = trim(dir_global)//adjustl(fn_global)
     call init_CT_lib( fng, ptop_wrf, xlon, xlat, znu, nx, ny, nz, day_start, hour_start, gincr )
     print*,'call to initialize global model'
  endif
  

  !   given the proper dimensions, allocate the boundary and tendency variables
  allocate( bcxs(ny,nz,nw) )
  allocate( bcxe(ny,nz,nw) )
  allocate( bcys(nx,nz,nw) )
  allocate( bcye(nx,nz,nw) )

  allocate( btxs(ny,nz,nw) )
  allocate( btxe(ny,nz,nw) )
  allocate( btys(nx,nz,nw) )
  allocate( btye(nx,nz,nw) )

  ! initialize tendencies as part of solution to problem
  ! encountered if only one time step in wrfbdy
  btxs = 0.
  btxe = 0.
  btys = 0.
  btye = 0.

  ! these only needed for chem_bc_opt 1
  if(chem_bc_opt == 1) then
     allocate( bcxso(ny,nz,nw) )
     allocate( bcxeo(ny,nz,nw) )
     allocate( bcyso(nx,nz,nw) )
     allocate( bcyeo(nx,nz,nw) )
  endif


  !---------------------------------------------------------------------------------------------------
  !     do the time stepping
  !---------------------------------------------------------------------------------------------------
 

  ! Variable 'it_bdy' represents the time increment in wrfbdy; 'it_glob' represents a time increment in 
  !   the global model 2-month consolidated files
  ! Note that time stepping logic is dependent on the time resolution of the global model source
  ! relative to the wrfbdy time step.

  do it_bdy = 1, ntime
     write(*,*) 'Time step number of wrfbdy : ',it_bdy 
      
        if(chem_bc_opt == 1) then

           if (it_bdy .eq. 1) then
             it_glob = (day_start - 1)*gincr+1
           else
             it_glob = (day_start - 1)*gincr+(it_bdy-1)*2
           endif
           if (hour_start == 12) then
             it_glob = it_glob + gincr/2
           endif
 
           write(*,*) 'Time increment of global model  : ',it_glob

           call global_interpolate4d( specname, bcxs, bcxe, bcys, &
                bcye, nx, ny, nz, nw, it_bdy, it_glob, dt_w, dt_g )
        endif
 
        !  use option 2 to avoid losing sink flux signals (subtract this value in postprocessing)          
        if(chem_bc_opt == 2) then 
           bcxs(:,:,:)=300.
           bcxe(:,:,:)=300.
           bcys(:,:,:)=300.
           bcye(:,:,:)=300.
        endif

        ! use option 3 only for tracers with positive definite flux       
        if(chem_bc_opt == 3) then 
           bcxs(:,:,:)=0.
           bcxe(:,:,:)=0.
           bcys(:,:,:)=0.
           bcye(:,:,:)=0.
        endif
      

        ! [PSU] tendencies for options 2 and 3 updated December 2015
        !       These were set at 1e-5 which is not nearly small enough,
        !       and introduced 0.036 ppm 'leaks' per tracer in waves on inflows on the walls,
        !       for an annual run, this amounts to approximately 0.1 ppm per tracer (option 2 and 3) in the interior.
        !       Tendencies changed to 0, which seems to be interpreted as approximately 1.e-15.

        if (it_bdy >=2) then

           ! Compute Tendencies
           if(chem_bc_opt == 1) then 
              btxs = ( bcxs - bcxso )/(float(dt_w)*3600.)
              btxe = ( bcxe - bcxeo )/(float(dt_w)*3600.)
              btys = ( bcys - bcyso )/(float(dt_w)*3600.)
              btye = ( bcye - bcyeo )/(float(dt_w)*3600.)
           endif

           if(chem_bc_opt == 2) then 
              btxs = 0.
              btxe = 0.
              btys = 0.
              btye = 0.
           endif

           if(chem_bc_opt == 3) then 
              btxs = 0.
              btxe = 0.
              btys = 0.
              btye = 0.
           endif

        endif

        ! save current bc's for next tendency computationn
        if (chem_bc_opt == 1) then
           bcxso = bcxs
           bcxeo = bcxe
           bcyso = bcys
           bcyeo = bcye
        end if
      
        write(*,*) 'Write the boundary conditions variables in wrfbdy file for time increment ', it_bdy
        call wrfchem_write4d( specname, bcxs, bcxe, bcys, bcye, btxs, btxe, btys, btye, it_bdy )
    
   
  end do

  !     exit from libs
  call exit_wrfchem_lib( )
  if(chem_bc_opt == 1) then
     call exit_global_lib( )
  endif

  !     deallocate memory space
  if( allocated( znu ) ) deallocate( znu )
  if( allocated( xlon ) ) deallocate( xlon )
  if( allocated( xlat ) ) deallocate( xlat )
  if( allocated( bcxs ) ) deallocate( bcxs )
  if( allocated( bcxe ) ) deallocate( bcxe )
  if( allocated( bcys ) ) deallocate( bcys )
  if( allocated( bcye ) ) deallocate( bcye )
  if( allocated( btxs ) ) deallocate( btxs )
  if( allocated( btxe ) ) deallocate( btxe )
  if( allocated( btys ) ) deallocate( btys )
  if( allocated( btye ) ) deallocate( btye )
  if( allocated( bcxso ) ) deallocate( bcxso )
  if( allocated( bcxeo ) ) deallocate( bcxeo )
  if( allocated( bcyso ) ) deallocate( bcyso )
  if( allocated( bcyeo ) ) deallocate( bcyeo )


  write(*,*) 'bc_wrfchem completed successfully'
end program main_bc_wrfchem_CT
