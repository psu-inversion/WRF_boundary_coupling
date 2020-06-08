
      module module_CT_lib

      use netcdf      
      implicit none

!     include files
      include 'netcdf.inc'

!     public procedures
      public  :: init_CT_lib, global_interpolate4d, exit_global_lib

!     private procedures
      private :: hybrid_profile_TM5, handle_error_global
      !private :: hybrid_profile_TM5  

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
      ! use of the following is highly model dependent (not generally used with CarbonTracker)
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

      subroutine init_CT_lib( global_fn, ptop_wrf, x2d, y2d, znu, nx, ny, nz, day_start, hour_start, gincr )

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
      real, allocatable :: dateCT(:,:)
      real, allocatable :: press_edges(:)                  ! used with hybrid_profile subroutine

      write(*,*) 'open the input netCDF file ',global_fn
      status = nf90_open( global_fn, nf90_nowrite, ncidg )
      if( status /= nf90_noerr ) call handle_error_global( status )

      write(*,*) 'get the longitude dimension length of Carbon Tracker'
      status = nf90_inq_dimid( ncidg, 'lon', dimid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, dimid, len=nlon )
      if( status /= nf90_noerr )  call handle_error_global( status )

      write(*,*) 'get the latitude dimension length of the Carbon Tracker'
      status = nf90_inq_dimid( ncidg, 'lat', dimid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, dimid, len=nlat )
      if( status /= nf90_noerr )  call handle_error_global( status )

!     get the vertical dimension length of the Carbon Tracker data
      status = nf90_inq_dimid( ncidg, 'level', dimid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, dimid, len=nlev )
      if( status /= nf90_noerr )  call handle_error_global( status )

      write(*,*) 'read longitudes'
      allocate( x1d(nlon), y1d(nlat) )

      status = nf90_inq_varid( ncidg, 'lon', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_get_var( ncidg, varid, x1d )
      if( status /= nf90_noerr )  call handle_error_global( status )
      ! if  the global model longitude vector is not -180 to +180, 
      !    then fix it here to match the wrf protocol

      status = nf90_inq_varid( ncidg, 'lat', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_get_var( ncidg, varid, y1d )
      if( status /= nf90_noerr )  call handle_error_global( status )

      write(*,*) 'read number of dates in the file'
      status = nf90_inq_dimid( ncidg, 'date', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, varid, len=ndate )
      if( status /= nf90_noerr ) call handle_error_global( status )

      write(*,*) 'read date format'
      status = nf90_inq_dimid( ncidg, 'calendar_components', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_inquire_dimension( ncidg, varid, len=ncalcomp )
      if( status /= nf90_noerr ) call handle_error_global( status )

      allocate ( dateCT(ncalcomp,ndate) )

      write(*,*) 'check time dimensions : ',ncalcomp, ndate

      status = nf90_inq_varid( ncidg, 'date_components', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      status = nf90_get_var( ncidg, varid, dateCT )
      if( status /= nf90_noerr )  call handle_error_global( status )

      ! dates are used strictly for reporting content of consolidated file  
      
      write(*,*) 'Dates in CT file - Start : ',dateCT(:,1)
      write(*,*) 'Dates in CT file - End : ',dateCT(:,ndate)

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

      write(*,*) 'x2d : ',x2d(1:10,23)  !debugging

      ! all domain
      do i = 1, nx
      do j = 1, ny 
         do n = 1, nlon
            if( x2d(i,j) < x1d(n) ) exit
         end do
         
         ix(0,i,j) = min( nlon-1, max(n-1,1) )
         ix(1,i,j) = ix(0,i,j) + 1
         ax(0,i,j) = ( x2d(i,j) - x1d(ix(0,i,j)) )/( x1d(ix(1,i,j)) - x1d(ix(0,i,j)) )
         ax(1,i,j) = 1.0 - ax(0,i,j)
         
         do n = 1, nlat
            if( y2d(i,j) < y1d(n) ) exit
         end do
         
         jy(0,i,j) = min( nlat-1, max(n-1,1) )
         jy(1,i,j) = jy(0,i,j) + 1
         by(0,i,j) = ( y2d(i,j) - y1d(jy(0,i,j)) )/( y1d(jy(1,i,j)) - y1d(jy(0,i,j)) )
         by(1,i,j) = 1.0 - by(0,i,j)
         
      end do
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

      status = nf90_inq_varid( ncidg, 'pressure', varid )
      if( status /= nf90_noerr )  call handle_error_global( status )

      nstt(1:4) = (/ 1, 1, 1, start_index /)   
      ncnt(1:4) = (/ nlon, nlat, 1 , 1/)
      status = nf90_get_var( ncidg, varid, ps_glob, nstt(1:4), ncnt(1:4) )
      if( status /= nf90_noerr ) call handle_error_global( status )

      !for clarification, ps_glob is CT surface pressure and
      !                   ps_wrf is CT surface pressure interpolated to the WRF grid 
      ! vertical level mapping is done with pressure columns in units of Pa
      ! both WRF and CT pressure units are in Pa
    
      
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

      write(*,*) 'vertical interpolation coefs'
      allocate( p_glob(nlev), p_wrf(nz) )
      p_glob(:)=0.
      p_wrf(:)=0.

      allocate( press_edges(nlev) )  !minus the topmost level 
      press_edges(:) = 0.


      !   At a minimum, only kz[0,i,j,k] is used in the interpolate subroutine
      !      NO horizontal interpolation of co2 field
      !   Computation of kz[0:1,i,j,k] and the corresponding czxs,czxe,czys,czye vertical
      !   interpolation coefficients done here so that they may be used in interpolate subroutine
      !   (if needed)

      !west
      i = 1
      do j = 1, ny

         ! The call to hybrid_profile_TM5 returns a column of pressures on the edges of 
         !    global model levels based on the interpolated surface pressure ps_wrf(i,j).
         ! The assignment of source global model levels to receiving wrf model levels
         !    is determined by the pressure edges between which the wrf level mid-point
         !    pressure falls.  This determination is done in log space with the respective
         !    pressure columns in units of Pa.

         call hybrid_profile_TM5( ps_wrf(i,j), nlev, press_edges )
         
         p_glob(:) = log( press_edges(:) )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf )
         
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
        ! more debugging displays
        write(*,*) 'interpolated surface pressure: ',ps_wrf(1,ny)
        write(*,*) 'pressure WRF : ',p_wrf(:)
        write(*,*) 'pressure CT : ',p_glob(:)
        write(*,*) 'kz levels : ',kz(0,1,ny,:)
 

      ! east
      i = nx
      do j = 1, ny

         call hybrid_profile_TM5( ps_wrf(i,j), nlev, press_edges )
         
         p_glob(:) = log( press_edges(:) )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf )
         
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

         call hybrid_profile_TM5( ps_wrf(i,j), nlev, press_edges )
         
         p_glob(:) = log( press_edges(:) )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf )
         
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

         call hybrid_profile_TM5( ps_wrf(i,j), nlev, press_edges )
         
         p_glob(:) = log( press_edges(:) )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf )

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
      deallocate( dateCT, x1d, y1d, ps_glob, ps_wrf, p_glob, p_wrf, press_edges )

      write(*,*) 'allocate memory space for reading'  
      allocate( globval(nlon,nlat,nlev), tmpval(nlon,nlat,nlev) )

      write(*,*) 'End of subroutine init_CT_lib'
      end subroutine init_CT_lib

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

!     output arguments
      real, intent(inout) :: wrfxs(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(inout) :: wrfxe(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(inout) :: wrfys(nx,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(inout) :: wrfye(nx,nz,nw)         ! wrfchem vmr(ppm)

!     local arguments
      integer :: status, co2_varid, new_itg, incrtoavg 
     
      character(20) :: globspn    ! variable name for co2 in the global file

      integer :: i, j, k, m

      ! not exactly sure why this is done (initialize to something/anything?)
      wrfxs(:,:,:) = 385.
      wrfxe(:,:,:) = 385.
      wrfys(:,:,:) = 385.
      wrfye(:,:,:) = 385.


!     read CarbonTracker CO2 mole fractions
!            'it_bdy' (wrfbdy time increment) and 'it_glob' (increment in global file) in main_bc_wrfchem
!             become known as 'itb' and 'itg' here in the global module



      globval(:,:,:) = 0.   
      globspn = 'co2'

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
           
            write(*,*) 'Read and interpolate Carbon Tracker CO2 : ', trim(globspn), new_itg
           
            nstt(1:4) = (/ 1, 1, 1, new_itg /)
            ncnt(1:4) = (/ nlon, nlat, nlev, 1 /)
            status = nf90_get_var( ncidg, co2_varid, tmpval(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf90_noerr ) call handle_error_global( status )
            
            globval(:,:,:) = globval(:,:,:) + tmpval(:,:,:)
   
            !debugging         
            !write(*,*) 'Some value just after reading the Carbon Tracker CO2 : ', globval(23:26,23,1)
        
      enddo
     
      globval(:,:,:) = globval(:,:,:)/float(incrtoavg)   !average three or six hours of global model
    

      write(*,*) 'Examples final :: ', globval(23,23,1),globval(3,4,5)

!     assign co2 by level in the wrf boundary grid cells
!     Note that only ix(0,i,j), jy(0,i,j), are used here horizontally, and kz(0,i,j,k) vertically
!     If considering using the vertical interpolation coefficients cz.., test to see if this is a good choice. 
!         (found not to be optimal for CarbonTracker)
!     If cz coefficients are used, uncomment the first 2 lines, and comment out the following line for
!            each of the cardinal directions
      
      ! west
      i= 1
      do k = 1, nz
         do j = 1, ny
            ! if using the cz...vertical interpolation coefficients, use these two lines:
            !wrfxs(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czxs(1,i,j,k) + &
            !                globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czxs(0,i,j,k)

            ! OR if not using the vertical interpolation coefficients (the usual case), use this line:
            wrfxs(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))
         enddo
      enddo
      

      ! east
      i = nx
      do k = 1, nz
         do j = 1, ny
            !wrfxe(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czxe(1,i,j,k) + &
            !                globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czxe(0,i,j,k)
            
            wrfxe(j,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))

         enddo
      enddo

      ! north
      j = ny
      do k = 1, nz
         do i = 1, nx
            !wrfye(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czye(1,i,j,k) + &
            !                globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czye(0,i,j,k)

            wrfye(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))
            
         enddo
      enddo

      ! south
      j = 1
      do k = 1, nz
         do i = 1, nx
            !wrfys(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czys(1,i,j,k) + &
            !                globval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czys(0,i,j,k) 

            wrfys(i,k,:)  = globval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))

         enddo
      enddo

      end subroutine global_interpolate4d

!---------------------------------------------------------------------------------------------------
!     hybrid_profile_TM5 - compute level edges of pressure column given surface pressure
!---------------------------------------------------------------------------------------------------

     subroutine hybrid_profile_TM5( psurf, nlevel, pcolumn_edges )

     ! Given a surface pressure (input) and the at and bt coefficients for the TM5 global model
     ! generate a pressure column on TM5 levels, from surface to model top.
     ! Return the pressures at nlevel boundary edges (omitting top edge).
     ! Use in the calling program stops well below the top edge.
     ! algorithm:  P = at + bt*psurf

     implicit none

     !input arguments:

     real, intent(in)                       :: psurf
     integer, intent(in)                    :: nlevel 
     real, dimension(nlevel),intent(inout)  :: pcolumn_edges    

     !local variables:

     integer                                :: kk
     integer, parameter                     :: boundary = 26
     real, parameter, dimension(boundary)   :: at = &
            (/0.0000000e+00, 7.3677430e+00, 2.1039389e+02, 8.5536176e+02, &
              2.0637798e+03, 3.8509133e+03, 6.1443149e+03, 8.8023564e+03, &
              1.1632759e+04, 1.4411124e+04, 1.6899469e+04, 1.8864750e+04, &
              2.0097402e+04, 2.0429863e+04, 1.9755109e+04, 1.8045184e+04, &
              1.5379806e+04, 1.2077446e+04, 8.7650537e+03, 6.0180195e+03, &
              3.9602915e+03, 1.6806403e+03, 7.1321808e+02, 2.9849579e+02, &
              9.5636963e+01, 0.0000000e+00/)    
     real, parameter, dimension(boundary)   :: bt = &
            (/1.0000000e+00, 9.9401945e-01, 9.7966272e-01, 9.5182151e-01, &
              9.0788388e-01, 8.4737492e-01, 7.7159661e-01, 6.8326861e-01, &
              5.8616841e-01, 4.8477158e-01, 3.8389215e-01, 2.8832296e-01, &
              2.0247590e-01, 1.3002251e-01, 7.3533830e-02, 3.4121160e-02, &
              1.1142910e-02, 1.8151600e-03, 7.5820000e-05, 0.0000000e+00, &
              0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, &
              0.0000000e+00, 0.0000000e+00/)    



     ! skipping the upper edge (which is at 0.0):
     do kk = 1, nlevel        
         pcolumn_edges(kk) = at(kk) + bt(kk) * psurf
     end do

     end subroutine hybrid_profile_TM5

!---------------------------------------------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------------------------------------------

      subroutine handle_error_global( status )

      implicit none

!     input arguments :
      integer, intent(in) :: status

!     print the error information from processing NETcdf file
      print*, nf90_strerror( status )

!     exit from the global lib
      call exit_global_lib( flag=1 )

      end subroutine handle_error_global

!---------------------------------------------------------------------------------------------------
!     exit from module_CT_lib 
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
        stop ' in module_CT_lib ...'
      else
        print*, 'successfully exit from module_CT_lib ...'
      endif

      end subroutine exit_global_lib


      end module module_CT_lib
