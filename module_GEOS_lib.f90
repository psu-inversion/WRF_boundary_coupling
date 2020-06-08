
      module module_GEOS_lib

      use netcdf      
      implicit none

!     include files
      include 'netcdf.inc'

!     public readprocedures
      public  :: init_GEOS_lib, global_interpolate4d, exit_global_lib

!     private procedures
      private :: hybrid_profile_GEOS, handle_error_global

!     public variables   check this: these may be used only in this module, so private?
      integer, public  :: nlon, nlat, nlev                    ! 3-D data dimensions of global model grid

!     variables for global model netCDF data
      integer, private :: ncidg                              ! netCDF file ID
      integer, private :: nstt(4), ncnt(4)                   ! start index and counter array
      integer, private :: start_index                        ! index to surface pressure to be interpolated
      integer, private :: ndate                              ! number of records in global consolidated file
      integer, private :: ncalcomp                           ! used to interpret/report the time span of the global file


!     variables used by interpolation, created in init subroutine and used in interpolate subroutine
!     [PSU]:  note that the only variables currently used are ix, jy, kz, ax, by,
!             and the cz.. series (czxs, czxe, czys, czye) is used only if warranted (test it!)

      integer, private, allocatable :: ix(:,:,:), jy(:,:,:), kz(:,:,:,:)  ! indices used in interpolation
      real, private, allocatable :: ax(:,:,:), by(:,:,:)                  ! horizontal weight coefficients
      ! use of the following is highly model dependent (used for CMS GEOS-Chem)
      real, private, allocatable :: czxs(:,:,:,:)     ! vertical weight coef. west
      real, private, allocatable :: czxe(:,:,:,:)     ! vertical weight coef. east
      real, private, allocatable :: czys(:,:,:,:)     ! vertical weight coef. south
      real, private, allocatable :: czye(:,:,:,:)     ! vertical weight coef. north

      ! used only in interpolate routine, but allocated during initialization

      real, private, allocatable :: globval(:,:,:), tmpval(:,:,:)

      contains

!---------------------------------------------------------------------------------------------------
!     initialize netcdf file
!---------------------------------------------------------------------------------------------------

      subroutine init_GEOS_lib( global_fn, ptop_wrf, x2d, y2d, znu, nx, ny, nz, day_start, hour_start, gincr )

      implicit none

!     input arguments
      integer, intent(in)      :: nx, ny, nz               !dimensions of wrf grid
      character(len=*), intent(in) :: global_fn 
      real, intent(in)         :: ptop_wrf                 !model top from wrf
      real, intent(in)         :: x2d(nx,ny), y2d(nx,ny)   !xlong and xlat from wrf
      real, intent(in)         :: znu(nz)                  !vertical discretization from wrf
      integer, intent(in)      :: day_start
      integer, intent(in)      :: hour_start
      integer, intent(in)      :: gincr         

!     local arguments
      integer :: status
      integer :: dimid, varid, co2_varid
      integer :: i, j, k, n
      real, allocatable :: x1d(:), y1d(:)                  ! longitude and latitude vectors from global model
      real, allocatable :: p_glob(:), p_wrf(:)             ! temporary use vertical interpolation vectors
      real, allocatable :: ps_glob(:,:), ps_wrf(:,:)       ! global model surface pressure and interpolated to wrf grid
      real, allocatable :: dateGEOS(:)
      real, allocatable :: press_edges(:)                  ! used with hybrid_profile subroutine


      write(*,*) 'open the input netCDF file ',global_fn
      status = nf90_open( global_fn, nf90_nowrite, ncidg )
      if( status /= nf90_noerr ) call handle_error_global( status )

      write(*,*) 'get the longitude dimension length of GEOS-CHEM'
      status = nf90_inq_dimid( ncidg, 'lon', dimid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, dimid, len=nlon )
      if( status /= nf90_noerr )  call handle_error_global( status )

      write(*,*) 'get the latitude dimension length of the GEOS-CHEM'
      status = nf90_inq_dimid( ncidg, 'lat', dimid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, dimid, len=nlat )
      if( status /= nf90_noerr )  call handle_error_global( status )

!     get the vertical dimension length of the GEOS-Chem data
      status = nf90_inq_dimid( ncidg, 'lev', dimid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, dimid, len=nlev )
      if( status /= nf90_noerr )  call handle_error_global( status )


      write(*,*) 'read longitudes and latitudes'
      allocate( x1d(nlon), y1d(nlat) )

      status = nf90_inq_varid( ncidg, 'lon', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_get_var( ncidg, varid, x1d )
      if( status /= nf90_noerr )  call handle_error_global( status )
!     CMS GEOS-Chem longitude vector
!            in input file, it is [180.0, 185.0,...,355.0, 0.0, 5.0,...,175.0]
!     realign it here to match WRF protocol of [-180.0,...,0.0,...180.0]
      where( x1d > 179.0 )
        x1d = x1d - 360.0
      endwhere

!     read latitudes
      status = nf90_inq_varid( ncidg, 'lat', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_get_var( ncidg, varid, y1d )
      if( status /= nf90_noerr )  call handle_error_global( status )

      write(*,*) 'read number of dates in the file'
      status = nf90_inq_dimid( ncidg, 'time', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, varid, len=ndate )
      if( status /= nf90_noerr ) call handle_error_global( status )

      ! time in CMS GEOS-Chem file is hours since 1985 1-1 0 UTC
      ! report starting and ending times, number of days here
      allocate( dateGEOS(ndate) )

      status = nf90_get_var( ncidg, varid, dateGEOS )
      if( status /= nf90_noerr ) call handle_error_global( status )
    
      write(*,*) 'Starting hour in GEOS-Chem file:  ',dateGEOS(1)
      write(*,*) 'Ending hour in GEOS-Chem file  :  ',dateGEOS(ndate)
      write(*,*) 'Number of days in GEOS-Chem file : ', ndate/gincr

      ! issue warning (and fail) here if WRF domain is not entirely
      ! contained within global model domain

      if ( ( minval(x2d) < x1d(1) ) .or. &
           ( maxval(x2d) > x1d(nlon) ) .or. & 
           ( minval(y2d) < y1d(1) ) .or. &
           ( maxval(y2d) > y1d(nlat) ) ) then
          write(*,*) ' ***DOMAIN MISMATCH*** '
          call exit_global_lib( flag=2 )
      endif


!     allocate memory space to store interpolation coefficients
!      ax, by, ix, jy, kz, and (sometimes) cz.. 

      allocate( ax(0:1,nx,ny), by(0:1,nx,ny) )
      allocate( ix(0:1,nx,ny), jy(0:1,nx,ny) )
      allocate( kz(0:1,nx,ny,nz) )
      allocate( czxs(0:1,nx,ny,nz), czxe(0:1,nx,ny,nz), czys(0:1,nx,ny,nz), czye(0:1,nx,ny,nz) )



!     determine horizontal interpolation coefficients

      write(*,*) 'x2d : ',x2d(1:10,23)  !debugging display

      ! all domain
      do i = 1, nx
      do j = 1, ny 
         do n = 1, nlon
            if( x2d(i,j) < x1d(n) ) exit
         enddo
         
         ix(0,i,j) = min( nlon-1, max(n-1,1) )
         ix(1,i,j) = ix(0,i,j) + 1
         ax(0,i,j) = ( x2d(i,j) - x1d(ix(0,i,j)) )/( x1d(ix(1,i,j)) - x1d(ix(0,i,j)) )
         ax(1,i,j) = 1.0 - ax(0,i,j)
         
         do n = 1, nlat
            if( y2d(i,j) < y1d(n) ) exit
         enddo
         
         jy(0,i,j) = min( nlat-1, max(n-1,1) )
         jy(1,i,j) = jy(0,i,j) + 1
         by(0,i,j) = ( y2d(i,j) - y1d(jy(0,i,j)) )/( y1d(jy(1,i,j)) - y1d(jy(0,i,j)) )
         by(1,i,j) = 1.0 - by(0,i,j)
         
      enddo
   end do

      ! either set start_index to 1 here (old protocol)
      ! or find the time in global model file corresponding to day and hour start in wrf:
      
      ! start_index = 1
      !   OR
      start_index = (day_start -1) * gincr + 1
      if (hour_start == 12) then
         start_index = start_index + gincr/2
      end if 

      write(*,*) 'using global surface pressure from time step ', start_index

      write(*,*) 'read surface pressure and interpolate ...'
      allocate( ps_glob(nlon,nlat), ps_wrf(nx,ny) )

      status = nf90_inq_varid( ncidg, 'PEDGE_S__PSURF', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      nstt(1:3) = (/ 1, 1, start_index /)  
      ncnt(1:3) = (/ nlon, nlat, 1/)
      status = nf90_get_var( ncidg, varid, ps_glob, nstt(1:3), ncnt(1:3) )
      if( status /= nf90_noerr ) call handle_error_global( status )

      ! for clarification, ps_glob is GEOS-Chem surface pressure and
      !                   ps_wrf is GEOS-Chem surface pressure interpolated to the WRF grid 
      ! vertical level mapping is done with pressure columns in units of Pa
      ! Convert GEOS-Chem pressure units from hPa to Pa after computing pressure profile  
      !ps_glob(:,:) = ps_glob(:,:)*1e2

      ! use this surface pressure map to create interpolated surface pressure map in wrf coordinates  
      !write(*,*) start_index,ps_moz(23,23)

      write(*,*) 'Test horizontal interpolation coefficients : ',ix(0,23,23),ax(0,23,23),by(0,23,23),jy(0,23,23)
      do j = 1, ny 
      do i = 1, nx 
        ps_wrf(i,j) = (ps_glob(ix(0,i,j),jy(0,i,j))*ax(1,i,j)*by(1,i,j) +   &
                       ps_glob(ix(0,i,j),jy(1,i,j))*ax(1,i,j)*by(0,i,j) +   &
                       ps_glob(ix(1,i,j),jy(0,i,j))*ax(0,i,j)*by(1,i,j) +   &
                       ps_glob(ix(1,i,j),jy(1,i,j))*ax(0,i,j)*by(0,i,j))
      enddo
      enddo


      write(*,*) 'Examples of interpolated surface pressure : ',ps_wrf(1:5,2)

      write(*,*) 'vertical interpolation coefficients'
      allocate( p_glob(nlev), p_wrf(nz) )
      p_glob(:)=0.
      p_wrf(:)=0.

      allocate( press_edges(nlev) )  !minus the topmost level (used by hybrid profile)  
      press_edges(:) = 0.


      !   At a minimum, only kz[0,i,j,k] is used in the interpolate subroutine
      !      NO horizontal interpolation of co2 field
      !   Computation of kz[0:1,i,j,k]; the corresponding czxs,czxe,czys,czye vertical
      !   interpolation coefficients done here so that they may be used in interpolate subroutine
      !   (if needed)

      
      ![PSU]: note that before September 2015, only kz[0,i,j,k] is used in the interpolate subroutine
      !   This has been revised to use kz[0:1,i,j,k] and the corresponding czxs,czxe,czys,czye vertical
      !   interpolation coefficients
      !west
      i = 1
      do j = 1, ny

         ! The call to hybrid_profile_GEOS returns a column of pressures on the edges of 
         !    global model levels based on the interpolated surface pressure ps_wrf(i,j).
         !    This is done with pressure units in hPa.
         ! The assignment of source global model levels to receiving wrf model levels
         !    is determined by the pressure edges between which the wrf level mid-point
         !    pressure falls.  This determination is done in log space with the respective
         !    pressure columns in units of Pa.

         call hybrid_profile_GEOS( ps_wrf(i,j), nlev, press_edges )
 
         p_glob(:) = log( press_edges(:)*1e2 )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j)*1e2 + (1.0-znu(:))*ptop_wrf )

         do k = 1, nz
            
            do n = 1, nlev 
               if( p_wrf(k) > p_glob(n) ) exit
            enddo
            
            kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
            kz(1,i,j,k) = kz(0,i,j,k) + 1
            czxs(0,i,j,k) = ( p_wrf(k) - p_glob(kz(0,i,j,k)) )/( p_glob(kz(1,i,j,k)) - p_glob(kz(0,i,j,k)) )
            czxs(1,i,j,k) = 1.0 - czxs(0,i,j,k)
 
         enddo

          
      end do
         ! some debugging displays
        write(*,*) 'interpolated surface pressure: ',ps_wrf(1,ny)
        write(*,*) 'vertical interpolation vector for WRF : ',p_wrf(:)
        write(*,*) 'vertical interpolation vector for CMS : ',p_glob(:)
        write(*,*) 'kz levels : ',kz(0,1,ny,:)  
 

      ! east
      i = nx
      do j = 1, ny
         
         call hybrid_profile_GEOS( ps_wrf(i,j), nlev, press_edges )
         
         p_glob(:) = log( press_edges(:)*1e2 )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j)*1e2 + (1.0-znu(:))*ptop_wrf )

         
         do k = 1, nz
            
            do n = 1, nlev 
               if( p_wrf(k) > p_glob(n) ) exit
            enddo

            kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
            kz(1,i,j,k) = kz(0,i,j,k) + 1
            czxe(0,i,j,k) = ( p_wrf(k) - p_glob(kz(0,i,j,k)) )/( p_glob(kz(1,i,j,k)) - p_glob(kz(0,i,j,k)) )
            czxe(1,i,j,k) = 1.0 - czxe(0,i,j,k)
            
         enddo
      enddo

      ! north
      j = ny
      do i = 1, nx
         
         call hybrid_profile_GEOS( ps_wrf(i,j), nlev, press_edges )
         
         p_glob(:) = log( press_edges(:)*1e2 )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j)*1e2 + (1.0-znu(:))*ptop_wrf )
         
         
         do k = 1, nz
            
            do n = 1, nlev 
               if( p_wrf(k) > p_glob(n) ) exit
            enddo
            
            kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
            kz(1,i,j,k) = kz(0,i,j,k) + 1
            czye(0,i,j,k) = ( p_wrf(k) - p_glob(kz(0,i,j,k)) )/( p_glob(kz(1,i,j,k)) - p_glob(kz(0,i,j,k)) )
            czye(1,i,j,k) = 1.0 - czye(0,i,j,k)
            
         enddo
         
      enddo
      
      ! south
      j = 1
      do i = 1, nx

         call hybrid_profile_GEOS( ps_wrf(i,j), nlev, press_edges )
         
         p_glob(:) = log( press_edges(:)*1e2 )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j)*1e2 + (1.0-znu(:))*ptop_wrf )


        do k = 1, nz

          do n = 1, nlev 
            if( p_wrf(k) > p_glob(n) ) exit
          enddo

          kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
          kz(1,i,j,k) = kz(0,i,j,k) + 1
          czys(0,i,j,k) = ( p_wrf(k) - p_glob(kz(0,i,j,k)) )/( p_glob(kz(1,i,j,k)) - p_glob(kz(0,i,j,k)) )
          czys(1,i,j,k) = 1.0 - czys(0,i,j,k)

        enddo

      enddo

      write(*,*) 'release memory space not used'
      deallocate( dateGEOS, x1d, y1d, ps_glob, ps_wrf, p_glob, p_wrf, press_edges )

      write(*,*) 'allocate memory space for reading'
      allocate( globval(nlon,nlat,nlev), tmpval(nlon,nlat,nlev) )

      write(*,*) 'End of subroutine init_CT_lib'
      end subroutine init_GEOS_lib

!---------------------------------------------------------------------------------------------------
!     interpolate four-dimensional field ... 
!---------------------------------------------------------------------------------------------------

      subroutine global_interpolate4d( wrfspn, wrfxs, wrfxe, wrfys, wrfye, nx, ny, nz, nw, itb, itg, dtw, dtg  )
      use netcdf
      implicit none

!     input arguments
      integer, intent(in)      :: nx, ny, nz, nw     ! dimensions
      integer, intent(in)      :: itb, itg           ! aka it_bdy and it_glob in main
      integer, intent(in)      :: dtw, dtg           ! aka dt_w and dt_g in main
      character(len=*), intent(in) :: wrfspn         ! wrfchem species' name



      ! output arguments
      real, intent(out) :: wrfxs(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out) :: wrfxe(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out) :: wrfys(nx,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out) :: wrfye(nx,nz,nw)         ! wrfchem vmr(ppm)

!     local arguments
      integer :: status, co2_varid, new_itg, incrtoavg 
     
      character(20) :: globspn    ! variable name for co2 in the global file

      integer :: i, j, k, m

      ! not exactly sure why this is done (initialize to something/anything?)
      wrfxs(:,:,:) = 385.
      wrfxe(:,:,:) = 385.
      wrfys(:,:,:) = 385.
      wrfye(:,:,:) = 385.


!     read GEOS-Chem CO2 mole fractions
!            'it_bdy' (wrfbdy time increment) and 'it_glob' (increment in global file) in main_bc_wrfchem
!             become known as 'itb' and 'itg' here in the global module

      globval(:,:,:) = 0.
      globspn = 'IJ_AVG_S__NOx'
   

      ! get the varid of co2 in the global model
      status = nf90_inq_varid( ncidg, globspn, co2_varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      if (itb .eq. 1) then    !first 6-hr wrfbdy from 3 hours (1 increment) of global model, thereafter it is 2 increments
        incrtoavg = (dtw/dtg)/2
      else
        incrtoavg = dtw/dtg
      endif

      do m = 1, incrtoavg
        
            tmpval(:,:,:)=0.
            new_itg = itg + m -1
           
            write(*,*) 'Read GEOS/CMS CO2 : ', trim(globspn), new_itg
           
            nstt(1:4) = (/ 1, 1, 1, new_itg /)
            ncnt(1:4) = (/ nlon, nlat, nlev, 1 /)
            status = nf90_get_var( ncidg, co2_varid, tmpval(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf90_noerr ) call handle_error_global( status )
            
            globval(:,:,:) = globval(:,:,:) + tmpval(:,:,:)*1e-3   ! change unit to ppm from ppb
          
            ! debugging display
            !write(*,*) 'Some values just after reading the CMS CO2 : ', globval(23:26,23,1)
       
      enddo

      globval(:,:,:) = globval(:,:,:)/float(incrtoavg)   !average three or six hours of global model


      write(*,*) 'Examples final :: ', globval(23,23,1),globval(3,4,5)

!     assign co2 by level in the wrf boundary grid cells
!     Note that only ix(0,i,j), jy(0,i,j), must be used here horizontally, and kz(0,i,j,k) vertically
!     If considering using the vertical interpolation coefficients cz.., test to see if this is a good choice. 

!     If cz coefficients are used, use the first 2 lines, and comment out the following line for
!            each of the cardinal directions

   
      ! west
      i= 1
      do k = 1, nz
         do j = 1, ny
           ! if using the cz...vertical interpolation coefficients, use these two lines:
            wrfxs(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czxs(1,i,j,k) + &
                            globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czxs(0,i,j,k)

           ! OR if not using the vertical interpolation coefficients, use this line:
           ! wrfxs(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))
         enddo
      enddo
      

      ! east
      i = nx
      do k = 1, nz
         do j = 1, ny
            wrfxe(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czxe(1,i,j,k) + &
                            globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czxe(0,i,j,k)

            !OR
            !wrfxe(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))
         enddo
      enddo

      ! north
      j = ny
      do k = 1, nz
         do i = 1, nx
            wrfye(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czye(1,i,j,k) + &
                            globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czye(0,i,j,k) 
            !OR
            !wrfye(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))
            
         enddo
      enddo

      ! south
      j = 1
      do k = 1, nz
         do i = 1, nx
            wrfys(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czys(1,i,j,k) + &
                            globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czys(0,i,j,k)

            !OR 
            !wrfys(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))

         enddo
      enddo

      end subroutine global_interpolate4d

!---------------------------------------------------------------------------------------------------
!     hybrid_profile_GEOS - compute level edges of pressure column given surface pressure
!---------------------------------------------------------------------------------------------------

     subroutine hybrid_profile_GEOS( psurf, nlevel, pcolumn_edges )

     ! Given a surface pressure (input) and the ap and bp coefficients for the GEOS5 vertical profile
     ! generate a pressure column on GEOS5 levels, from surface to model top.
     ! Return the pressures at nlevel boundary edges (omitting top edge).
     ! Use in the calling program stops well below the top edge.
     ! algorithm:  pedge(k) = ap(k) + bp(k) * psurf(i,j)

     ! The algorithm is as supplied with the README with the CMS GEOS-Chem CO2 files and has been
     ! verified with the GEOS-Chem Wiki at
     !  //wiki.seas.harvard.edu/geos-chem/GEOS-Chem vertical_grids/
     ! This is the 47-layer reduced vertical grid for GMAO GEOS-5  (and MERRA)

     implicit none

     !input arguments:

     real, intent(in)                       :: psurf
     integer, intent(in)                    :: nlevel 
     real, dimension(nlevel),intent(inout)  :: pcolumn_edges    

     !local variables:

     integer                                :: kk
     integer, parameter                     :: boundary = 48
     real, parameter, dimension(boundary)   :: ap = &
              (/0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, &
                1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01, &
                4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01, &
                7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, &
                1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02, &
                1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02, &
                2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, &
                2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02, &
                1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01, &
                7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01, &
                1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00, &
                6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02/)
               
     real, parameter, dimension(boundary)   :: bp = &
              (/1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, &
                9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01, &
                8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01, &
                7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, &
                6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01, &
                4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01, &
                2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, &
                6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09, &
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, &
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, &
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, &
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00/) 
 

     ! skipping the upper edge (which is at 0.01 hPa):
     do kk = 1, nlevel        
         pcolumn_edges(kk) = ap(kk) + bp(kk) * psurf
     end do

     end subroutine hybrid_profile_GEOS

!---------------------------------------------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------------------------------------------

      subroutine handle_error_global( status )

      implicit none

!     input arguments :
      integer, intent(in) :: status

!     print the error information from processing NETcdf file
      print*, nf90_strerror( status )

!     exit from the bconLib
      call exit_global_lib( flag=1 )

      end subroutine handle_error_global

!---------------------------------------------------------------------------------------------------
!     exit from module_GEOS_lib 
!---------------------------------------------------------------------------------------------------

      subroutine exit_global_lib( flag )

!     input arguments
      integer, optional, intent(in) :: flag

!     local arguments
      integer :: status

!     release memory space
      if( allocated( ix) ) deallocate( ix )
      if( allocated( jy) ) deallocate( jy )
      if( allocated( kz) ) deallocate( kz )
      if( allocated( ax) ) deallocate( ax )
      if( allocated( by) ) deallocate( by )

      if( allocated(czxs) ) deallocate( czxs )
      if( allocated(czxe) ) deallocate( czxe )
      if( allocated(czys) ) deallocate( czys )
      if( allocated(czye) ) deallocate( czye )

      if( allocated( globval) ) deallocate( globval )
      if( allocated( tmpval) ) deallocate( tmpval )

!     close netCDF file
      if( ncidg /= 0 ) status = nf90_close( ncidg )

!     output information
      if( present(flag) ) then
        select case( flag )
          case( 1 ); print*, 'fail to process netCDF file...'
          case( 2 ); print*, 'fail for domain mismatch...'
          case default; print*, 'unknown error(s) occurred ...'
        endselect
        stop ' in module_GEOS_lib ...'
      else
        print*, 'successfully exit from module_GEOS_lib ...'
      endif

      end subroutine exit_global_lib


      end module module_GEOS_lib
