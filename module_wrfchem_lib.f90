
      module module_wrfchem_lib

      use netcdf
      implicit none

!     include files
      include 'netcdf.inc'

!     public  procedures
      public  :: init_wrfchem_lib, exit_wrfchem_lib
      public  :: wrfchem_readscalar, wrfchem_read2d,  wrfchem_read3d
      public  :: wrfchem_write4d

!     private procedures
      private :: handle_error, define_variable

!     variables defined here and used elsewhere (day_start and hour_start)
      integer, public  :: nx, ny, nz, nw, ntime      ! data dimension lengths
      integer, public  :: year_start, month_start, day_start
      integer, public  :: hour_start, minute_start, second_start

      real, public    :: ptop_wrf
      real, public, allocatable :: znu(:)
      real, public, allocatable :: xlon(:,:), xlat(:,:)


!     variables for reading netCDF data in this module 
      integer          :: ncid, ncidi          ! netCDF IDs for wrfbdy and wrfinput

      integer, private :: nstt(4), ncnt(4)     ! start index and counter array
      integer, private :: nscalar, nstring     ! data dimension length


      contains

!---------------------------------------------------------------------------------------------------
!     initialize netcdf file
!---------------------------------------------------------------------------------------------------

      subroutine init_wrfchem_lib( wrfchem_bdy_fn, wrfchem_input_fn, specname, ntime, ptop_wrf, znu, xlon, xlat )

      implicit none

!     input arguments
      integer, intent(out)              ::  ntime
      character(len=*), intent(in)      ::  wrfchem_bdy_fn, wrfchem_input_fn
      character(len=*), intent(in)      ::  specname
      real, intent(out)                 ::  ptop_wrf
      real, allocatable, intent(out)    ::  znu(:)
      real, allocatable, intent(out)    ::  xlon(:,:), xlat(:,:)


!     local arguments   NOTE: any changes here?  ;like character len's?
      integer :: status
      integer :: xdimid,  ydimid, zdimid, wdimid, tdimid, sdimid
      integer :: ndims, dimids(4)
      integer :: vid
      integer :: ip, i, n
      character(len=100) :: dtstrings, dtstringe  !, dtstringf   ! used to parse wrf date(s)
      ! these are used only if need to define variables in wrfbdy file
      character(len=100)  :: spec1
      character(len=3)   :: order

      write(*,*) 'init_wrfchem_lib: specname ',specname  

      write(*,*) 'open the bdy NETCDF file : ',wrfchem_bdy_fn
      status = nf90_open( wrfchem_bdy_fn, nf90_write, ncid )
      if( status /= nf90_noerr ) call handle_error( status )

      write(*,*) 'open the input NETCDF file : ', wrfchem_input_fn
      status = nf90_open( wrfchem_input_fn, nf90_nowrite, ncidi )
      if( status /= nf90_noerr ) call handle_error( status )

      write(*,*) 'get spatial dimension lengths from wrfbdy'
      status = nf90_inq_dimid( ncid, 'west_east', xdimid )
      if( status /= nf90_noerr ) call handle_error( status )

      status = nf90_inquire_dimension( ncid, xdimid, len=nx )
      if( status /= nf90_noerr ) call handle_error( status )

      status = nf90_inq_dimid( ncid, 'south_north', ydimid )
      if( status /= nf90_noerr ) call handle_error( status )

      status = nf90_inquire_dimension( ncid, ydimid, len=ny )
      if( status /= nf90_noerr ) call handle_error( status )

      status = nf90_inq_dimid( ncid, 'bottom_top', zdimid )
      if( status /= nf90_noerr ) call handle_error( status )

      status = nf90_inquire_dimension( ncid, zdimid, len=nz )
      if( status /= nf90_noerr ) call handle_error( status )

      status = nf90_inq_dimid( ncid, 'bdy_width', wdimid )
      if( status /= nf90_noerr ) call handle_error( status )

      status = nf90_inquire_dimension( ncid, wdimid, len=nw )
      if( status /= nf90_noerr ) call handle_error( status )

      write(*,*) 'get the time dimension length'
      status = nf90_inq_dimid( ncid, 'Time', tdimid )
      if( status /= nf90_noerr )  call handle_error( status )

      status = nf90_inquire_dimension( ncid, tdimid, len=ntime )
      if( status /= nf90_noerr ) call handle_error( status )

      write(*,*) 'get the string dimension length'
      status = nf90_inq_dimid( ncid, 'DateStrLen', sdimid )
      if( status /= nf90_noerr )  call handle_error( status )

      status = nf90_inquire_dimension( ncid, sdimid, len=nstring )
      if( status /= nf90_noerr ) call handle_error( status )

      write(*,*) 'enter redefine mode'
      status = nf90_redef( ncid )
      if( status /= nf90_noerr ) call handle_error( status )

      ! Note that these variables should already exist if wrfbdy created by WRF-Chem

      !  define variables if they do not already exist      
      ndims = 4

      !NOTE: have not completely updated subroutine define_variable for netCDF 90
      ! 'order' corrected to agree with wrfbdy created by WRF-Chem
      !  But note that descriptions and units are kinda generic out of WRF-Chem
      
      ! west-east
      dimids(1:4) = (/ ydimid, zdimid, wdimid, tdimid /)
      ! west
      do i = 1, 2
            if(i==1) spec1 = trim(specname)//'_BXS'
            if(i==2) spec1 = trim(specname)//'_BTXS'
            status = nf90_inq_varid( ncid, trim(spec1), vid )
            if( status /= nf90_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               !order = 'YSZ'
               order = 'XSZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            endif
      end do

      ! east
      do i = 1, 2 
         if(i==1) spec1 = trim(specname)//'_BXE'
         if(i==2) spec1 = trim(specname)//'_BTXE'
            status = nf90_inq_varid( ncid, trim(spec1), vid )
            if( status /= nf90_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               !order = 'YEZ'
               order = 'XEZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            end if
      end do
      
      ! north-south
      dimids(1:4) = (/ xdimid, zdimid, wdimid, tdimid /)
      ! north
      do i = 1, 2 
         if(i==1) spec1 = trim(specname)//'_BYS'
         if(i==2) spec1 = trim(specname)//'_BTYS'
            status = nf90_inq_varid( ncid, trim(spec1), vid )
            if( status /= nf90_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               !order = 'XEZ'
               order = 'YSZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            end if
      end do

      ! south
      do i = 1, 2 
         if(i==1) spec1 = trim(specname)//'_BYE'
         if(i==2) spec1 = trim(specname)//'_BTYE'
            status = nf90_inq_varid( ncid, trim(spec1), vid )
            if( status /= nf90_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               !order = 'XSZ'
               order = 'YEZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            end if
      end do
 

!     end define mode
      status = nf90_enddef( ncid )
      if( status /= nf90_noerr ) call handle_error( status )

      ! NOTE: only need first? do not need second if hard code dt_w, do not need last if know ntime
      !   changed for reporting, read start and end time stamps and simply report for documentation
      write(*,*) 'read start and end date time stamps' ! and second time stamp'
      status = nf90_inq_varid( ncid, 'Times', vid )
      if( status /= nf90_noerr )  call handle_error( status )

      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ nstring, 1 /)
      status = nf90_get_var( ncid, vid, dtstrings, nstt(1:2), ncnt(1:2) )

      nstt(1:2) = (/ 1, ntime /)
      status = nf90_get_var( ncid, vid, dtstringe, nstt(1:2), ncnt(1:2) )

      !nstt(1:2) = (/ 1, 2 /)
      !status = nf90_get_var( ncid, vid, dtstringf, nstt(1:2), ncnt(1:2) )

     ! debugging:
      write(*,*) 'ntime ', ntime
      write(*,*) 'wrfbdy starting time stamp ', dtstrings
      write(*,*) 'wrfbdy ending time stamp   ', dtstringe

     ! parse the time variable string to determine start date and time
     ! write(*,*) 'get the location of _ '
      ip = index( dtstrings, '_' )

      !write(*,*) 'delete -, _ and : '
      do n = 1, nstring 
        if( (dtstrings(n:n)=='-') .or. (dtstrings(n:n)=='_') .or.  &
             dtstrings(n:n)==':' ) then
          dtstrings(n:n) = ' '
        endif
      !  if( (dtstringe(n:n)=='-') .or. (dtstringe(n:n)=='_') .or.  &
      !       dtstringe(n:n)==':' ) then
      !    dtstringe(n:n) = ' '
      !  endif
        !if( (dtstringf(n:n)=='-') .or. (dtstringf(n:n)=='_') .or.  &
        !     dtstringf(n:n)==':' ) then
        !  dtstringf(n:n) = ' '
        !endif
      enddo

      !write(*,*) 'read year, month and day'
      ! day_start and hour_start used elsewhere....
      read( dtstrings(:ip), * ) year_start, month_start, day_start
      read( dtstrings(ip:nstring), * ) hour_start, minute_start, second_start
      !read( dtstringe(:ip), * ) year_end, month_end, day_end
      !read( dtstringe(ip:nstring), * ) hour_end, minute_end, second_end
      !read( dtstringf(:ip), * ) year_first, month_first, day_first
      !read( dtstringf(ip:nstring), * ) hour_first, minute_first, second_first


    ! These following are only needed for chem_bc_opt = 1, but were being accessed
    !  unconditionally before,  no harm done...

    !     read top ref. pressure from wrf input data
    call wrfchem_readscalar( 'P_TOP', ptop_wrf )
    print*, 'read  top ref. pressure from wrf input data'

    !     read eta values on half (mass) levels from wrfchem
    allocate( znu(nz) )
    call wrfchem_read2d( 'ZNU', znu )
    print*, 'read eta values on half (mass) levels',znu(1)

    !     read longitudes and latitudes from wrfchem input data
    allocate( xlon(nx,ny), xlat(nx,ny) )
    call wrfchem_read3d( 'XLONG', xlon )
    call wrfchem_read3d( 'XLAT',  xlat )
    print*,'read longitudes and latitudes from wrfchem input data'
 
      
      write(*,*) 'succesfully init module_wrfchem_lib ...'

      end subroutine init_wrfchem_lib

!---------------------------------------------------------------------------------------------------
!     define variable
!     1. Define variable according input arguments
!     2. Define the attributes of the variable
!---------------------------------------------------------------------------------------------------

     !NOTE: this subroutine has not (yet) been updated to netCDF/Fortran90 conventions

      subroutine define_variable( vname, nf_type, ndims, dimids, order )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname
      integer, intent(in)           :: nf_type
      integer, intent(in)           :: ndims
      integer, intent(in)           :: dimids(ndims)
      character (len=*), intent(in) :: order

!     local arguments
      integer :: status
      character(len=20) :: string 
      integer :: vid

!     define variables
      status = nf_def_var( ncid, vname, nf_type, ndims, dimids, vid )
      if( status /= nf_noerr )  call handle_error( status )

!     define attributes
      status = nf_put_att_int( ncid, vid, 'FieldType', nf_int, 1, 104 )
      if( status /= nf_noerr )  call handle_error( status )

      string = trim(order) 
      status = nf_put_att_text( ncid, vid, 'MemoryOrder', len_trim(string), trim(string) )
      if( status /= nf_noerr )  call handle_error( status )

      string = '-' 
      status = nf_put_att_text( ncid, vid, 'description', len_trim(string), trim(string) ) 
      if( status /= nf_noerr )  call handle_error( status )

      string = '-'
      status = nf_put_att_text( ncid, vid, 'units', len_trim(string), trim(string) )
      if( status /= nf_noerr )  call handle_error( status )

      string = ' ' 
      status = nf_put_att_text( ncid, vid, 'stagger', len_trim(string), trim(string) )
      if( status /= nf_noerr )  call handle_error( status )

      end subroutine define_variable

!---------------------------------------------------------------------------------------------------
!     read scalar data 
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_readscalar( vname, val )   ! read scalar variable 

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname

!     output arugments
      real, intent(out)        :: val

!     local arguments
      integer :: status
      integer :: varid

!     get variable ID, if no such a variable, return
      status = nf90_inq_varid( ncidi, vname, varid )
      if( status /= nf90_noerr ) call handle_error( status )
      write(*,*) ' wrfchem_readscalar: get ',vname,varid

!     read value
      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ 1, 1 /)
      status = nf90_get_var( ncidi, varid, val)
      if( status /= nf90_noerr ) call handle_error( status )

      write(*,*)' wrfchem_readscalar: value ',val,status

      end subroutine wrfchem_readscalar


!---------------------------------------------------------------------------------------------------
!     read one-dimensional data
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_read2d( vname, val )   ! read 2d variable( a slice of data )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname

!     output arugments
      real, intent(out)        :: val(nz)

!     local arguments
      integer :: status
      integer :: varid

!     get variable ID, if no such a variable, return
      status = nf90_inq_varid( ncidi, vname, varid )
      if( status /= nf90_noerr ) call handle_error( status )

!     read value
      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ nz, 1 /)
      status = nf90_get_var( ncidi, varid, val(:), nstt(1:2), ncnt(1:2) )
      if( status /= nf90_noerr ) call handle_error( status )

      end subroutine wrfchem_read2d

!---------------------------------------------------------------------------------------------------
!     read three-dimensional data
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_read3d( vname, val )   ! read 3d variable( a 2D array of data )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname

!     output arugments
      real, intent(out)        :: val(nx,ny)

!     local arguments
      integer :: status
      integer :: varid

!     get variable ID, if no such a variable, return
      status = nf90_inq_varid( ncidi, vname, varid )
      if( status /= nf90_noerr ) call handle_error( status )

!     read value
      nstt(1:3) = (/ 1, 1, 1 /)
      ncnt(1:3) = (/ nx, ny, 1 /)
      status = nf90_get_var( ncidi, varid, val(:,:), nstt(1:3), ncnt(1:3)  )
      if( status /= nf90_noerr ) call handle_error( status )

      end subroutine wrfchem_read3d

!---------------------------------------------------------------------------------------------------
!     write four-dimensional data
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_write4d( vname, valxs, valxe, valys, valye, vatxs, vatxe, vatys, vatye,it )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname
      real, intent(in)         :: valxs(ny,nz,nw)
      real, intent(in)         :: valxe(ny,nz,nw)
      real, intent(in)         :: valys(nx,nz,nw)
      real, intent(in)         :: valye(nx,nz,nw)
      real, intent(in)         :: vatxs(ny,nz,nw)
      real, intent(in)         :: vatxe(ny,nz,nw)
      real, intent(in)         :: vatys(nx,nz,nw)
      real, intent(in)         :: vatye(nx,nz,nw)
      integer, intent(in)      :: it

!     output arugments

!     local arguments
      integer :: status
      integer :: varid
      integer :: i
      character(len=100) :: bcname
      character(len=5),DIMENSION(8)  :: varext

      varext = (/'_BXS ', '_BTXS', '_BXE ', '_BTXE', '_BYS ', '_BTYS',&
                    '_BYE ', '_BTYE'/)

      do i = 1, 8
         bcname = trim(vname)//trim(varext(i))
         !write(*,*) 'get variable id in the file : ',bcname
         status = nf90_inq_varid( ncid, trim(bcname), varid )
         if( status /= nf90_noerr )  call handle_error( status )
         
         !write(*,*) 'set start and count arrays'
         nstt(1:4) = (/ 1, 1, 1, it /)
        
         if(i <= 4) ncnt(1:4) = (/ ny, nz, nw, 1 /)
         if(i >  4) ncnt(1:4) = (/ nx, nz, nw, 1 /)
         
         !write(*,*) 'write BC values'
         ! west
         if ( i == 1 ) then
            status = nf90_put_var( ncid, varid, valxs(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf90_noerr ) call handle_error( status )     
            ! east       
         elseif ( i == 3 ) then
            status = nf90_put_var( ncid, varid, valxe(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf90_noerr ) call handle_error( status )  
            ! north
         elseif ( i == 5 ) then
            status = nf90_put_var( ncid, varid, valys(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf90_noerr ) call handle_error( status ) 
            ! south
         elseif ( i == 7 ) then
            status = nf90_put_var( ncid, varid, valye(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf90_noerr ) call handle_error( status )
         else
            if ( it >= 2 ) then
               nstt(1:4) = (/ 1, 1, 1, it-1 /)
               ! Tendencies 
               ! west  
               if ( i == 2 ) then
                  status = nf90_put_var( ncid, varid, vatxs(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf90_noerr ) call handle_error( status )  
                  ! east        
               elseif ( i == 4 ) then
                  status = nf90_put_var( ncid, varid, vatxe(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf90_noerr ) call handle_error( status ) 
                  ! north
               elseif ( i == 6 ) then
                  status = nf90_put_var( ncid, varid, vatys(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf90_noerr ) call handle_error( status )  
                  ! south
               elseif ( i == 8 ) then
                  status = nf90_put_var( ncid, varid, vatye(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf90_noerr ) call handle_error( status )
               end if
            end if
         end if
         
      end do

      ! Tendencies above are written out offset from time step (for example, tendencies for
      ! time step 1 are written out at same time as boundaries for time step 2).
      ! To solve the problem of initializing tendencies for time step 1
      ! if there is only 1 time step in wrfbdy:
      ! Write the tendencies for time step 1 with the initial values supplied from
      ! the main code unconditionally at time step 1.  For applications with more than
      ! one time step in wrfbdy (common use), these will be overwritten at time step 2.

      if ( it == 1 ) then
         nstt(1:4) = (/ 1, 1, 1, it /)
         do i = 2, 8, 2
            !write(*,*) 'set start and count arrays'        
            if(i <= 4) ncnt(1:4) = (/ ny, nz, nw, 1 /)
            if(i >  4) ncnt(1:4) = (/ nx, nz, nw, 1 /)

            bcname = trim(vname)//trim(varext(i))
            !write(*,*) 'get variable id in the file : ',bcname
            status = nf90_inq_varid( ncid, trim(bcname), varid )
            if( status /= nf90_noerr )  call handle_error( status )

            ! Tendencies 
            ! west  
            if ( i == 2 ) then
               status = nf90_put_var( ncid, varid, vatxs(:,:,:), nstt(:), ncnt(:) )
               if( status /= nf90_noerr ) call handle_error( status )  
               ! east        
            elseif ( i == 4 ) then
               status = nf90_put_var( ncid, varid, vatxe(:,:,:), nstt(:), ncnt(:) )
               if( status /= nf90_noerr ) call handle_error( status ) 
               ! north
            elseif ( i == 6 ) then
               status = nf90_put_var( ncid, varid, vatys(:,:,:), nstt(:), ncnt(:) )
               if( status /= nf90_noerr ) call handle_error( status )  
               ! south
            elseif ( i == 8 ) then
               status = nf90_put_var( ncid, varid, vatye(:,:,:), nstt(:), ncnt(:) )
               if( status /= nf90_noerr ) call handle_error( status )
            end if
         end do
       end if


      end subroutine wrfchem_write4d 


!---------------------------------------------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------------------------------------------

      subroutine handle_error( status )

      implicit none

!     input arguments :
      integer, intent(in) :: status

!     print the error information from processing netCDF file
      print*, nf90_strerror( status )

!     exit from the wrfchem_lib
      call exit_wrfchem_lib( flag=1 )

      end subroutine handle_error

!---------------------------------------------------------------------------------------------------
!     exit from wrfchem_lib
!---------------------------------------------------------------------------------------------------

      subroutine exit_wrfchem_lib( flag )

!     input arguments
      integer, optional, intent(in) :: flag

!     local arguments
      integer :: status

!     close netCDF files
      if( ncid /= 0 ) status = nf90_close( ncid )
      if( ncidi /= 0 ) status = nf90_close( ncidi )

!     output information    NOTE:  I think this flag is set only in handle_error above?
      if( present(flag) ) then
        select case( flag )
          case( 1 ); print*, 'module_wrfchem_lib: fail to process netCDF file...'
          case( 2 ); print*, 'no such a species to save ...' 
          case default; print*, 'unknown error(s) occurred ...'
        endselect
        stop ' in module_wrfchem_lib ...'
      else
        print*, 'successfully exit from module_wrfchem_lib ...'
      endif 

      end subroutine exit_wrfchem_lib

      end module module_wrfchem_lib

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
