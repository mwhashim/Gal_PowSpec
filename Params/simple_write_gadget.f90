!===================================================
!
! Convert a single ascii-file to gadget-format
! 
! Set all the parameters in the subroutine
! gadget_from_ascii, compile and run!
!
!===================================================

! The header-values in a gadget-file
module gadget_header_params
  real(kind=8) :: massarr(0:5),a,redshift
  real(kind=8) :: omega0,omegaL,hubble,nall(0:5)
  real(kind=8) :: xLbox
  integer :: flag_sfr,flag_feedback,flag_cooling,numfiles
  integer :: nparttotal
end module gadget_header_params

module write_gadget_files
contains
  
  ! Write a gadget-file - assuming header-values are set in params
  ! and that pos,vel,mass,id are given
  subroutine write_gadget(filein,impi,jstart,jend,pos,vel,mass,id)
    use gadget_header_params
    implicit none
    character(len=*)    :: filein
    character(len=1024) :: filempi
    character(len=6)    :: string
    real(kind=4), dimension(:,:) :: pos
    real(kind=4), dimension(:,:) :: vel 
    real(kind=4), dimension(:) :: mass 
    real(kind=8) :: massarr_in(0:5),a_in,redshift_in
    real(kind=8) :: omega0_in,omegaL_in,hubble_in
    real(kind=8) :: xLbox_in
    integer, dimension(:) :: id
    integer :: npart_in(0:5), nall_in(0:5)
    integer :: lin, jstart, jend, i, j,impi
    integer :: flag_sfr_in,flag_feedback_in,flag_cooling_in,numfiles_in
    integer :: unused_in(64-6-12-2-2-1-1-6-1-1-2-2-2-2)

    ! Verbose
    write(*,*) 'Writing file i = ',impi

    ! Make output filename
    call convtoasc(impi,string)
    if (impi < 10) then
      filempi=trim(filein)//string(6:6)
    elseif (impi < 100) then
      filempi=trim(filein)//string(5:6)
    elseif (impi < 1000) then
      filempi=trim(filein)//string(4:6)
    else
      filempi=trim(filein)//string(3:6)
    endif

    ! Set header-values
    npart_in(0:5)    = 0
    npart_in(1)      = jend-jstart+1
    massarr_in       = massarr
    a_in             = a
    redshift_in      = redshift
    nall_in          = nall
    flag_sfr_in      = flag_sfr
    flag_feedback_in = flag_feedback
    nall_in          = nall
    flag_cooling_in  = flag_cooling
    numfiles_in      = numfiles
    xLbox_in         = xLbox
    omega0_in        = omega0
    omegaL_in        = omegaL
    hubble_in        = hubble

    ! Write header
    open(unit=lin,file=filempi,form='unformatted')
    write(lin) (npart_in(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
 &             redshift_in, flag_sfr_in, flag_feedback_in,  &
 &             (nall_in(i),i=0,5), flag_cooling_in, numfiles_in, &
 &             xLbox_in,omega0_in,omegaL_in,hubble_in,  &
 &             unused_in
      
    ! Write positons
    write(lin) ((pos(i,j),i=1,3),j=jstart,jend)
    
    ! Write velocity
    write(lin) ((vel(i,j),i=1,3),j=jstart,jend)
    
    ! Write id
    write(lin) (id(j),j=jstart,jend)

    ! If massarr(1) is set to 0 it means the masses are specified in
    ! the mass-array so write the mass array
    if(massarr_in(1) .eq. 0.0) then
      write(*,*) 'Writing masses'
      write(lin) (mass(j),j=jstart,jend)
    endif

  end subroutine write_gadget

  ! To convert an integer smaller than 999999 to a characters string
  subroutine convtoasc(number,sstring)

    implicit none
    integer :: number, istring, num, nums10, i
    character(len=6) :: sstring
    character(len=10), parameter :: nstring='0123456789'

    num=1000000
    nums10=num/10
    do i=1,6
       istring=1+mod(number,num)/nums10
       sstring(i:i)=nstring(istring:istring)
       num=num/10
       nums10=nums10/10
    enddo
  end subroutine convtoasc

  ! Takes an ascii-file and output a gadget file
  subroutine gadget_from_ascii()
    use gadget_header_params
    implicit none
    character(100) :: filein, fileout
    integer :: npart, lin, i, j
    real(kind=4), dimension(:,:), allocatable :: pos 
    real(kind=4), dimension(:,:), allocatable :: vel 
    real(kind=4), dimension(1:8) :: temp
    real(kind=4), dimension(:), allocatable :: mass
    integer, dimension(:), allocatable :: id
    real(kind=4) posunits
    real(kind=4) :: maxp, minp
    logical :: use_units_kpch   = .false.
    logical :: all_masses_equal = .true.

    ! In case we want pos/boxsize in the output file in units of kpc/h 
    if(use_units_kpch) then
      posunits = 1000.0
    else
      posunits = 1.0
    end if

    !=====================================================
    ! Set parameters
    !=====================================================

    ! Ascii-file with particle data
    filein = "tmpcat_761_laigle.dat"
    write(*,*) 'Opening ascii-file: ', trim(filein)
    
    ! Gadgetprefix for outputfile
    fileout = "gadget."

    ! Total number of particles
    npart        = 126361
    
    ! Set cosmological parameters
    a        = 1.0
    redshift = 0.0
    omega0   = 0.3
    omegaL   = 0.7
    hubble   = 0.704
    xLbox    = 100.0 * posunits

    ! Do not change these
    nall(0:5)    = 0
    massarr(0:5) = 0.0
    numfiles     = 1
    nparttotal   = npart
    nall(1)      = npart

    ! Mass of particles
    if(all_masses_equal) then
      massarr(1) = 3.0*omega0/dble(npart) * (xLbox/2998.0/posunits)**3 * 2.49265d11
      write(*,*) 'Setting all masses to be equal ', massarr(1)
      write(*,*) 'We will ignore the masses read from file'
    else
      massarr(1)   = 0.0
      write(*,*) 'We will read masses from file'
    end if

    ! Set flags (do not change these)
    flag_sfr      = 0
    flag_feedback = 0
    flag_cooling  = 0
 
    !========================================
    ! Read velocity, position and mass from file
    !========================================

    ! Allocate memory
    allocate(pos(3,npart))
    allocate(vel(3,npart))
    allocate(mass(npart))
    allocate(id(npart))

    pos(1:3,1:npart) = 0.0
    vel(1:3,1:npart) = 0.0
    mass(1:npart)    = 0.0
    id(1:npart)      = 0
   
    ! Open ascii-file and read pos, vel and mass
    open(unit=lin,file=filein,status='old',action='read')
    do i=1,npart
      ! Read: x, y, z, log mass, vx , vy, vz with x in Mpc, mass in Msun and v in km/s
      read(lin,*) temp(1:7)
     
      ! Set pos in units of Mpc/h (or kpc/h)
      pos(1:3,i) = (temp(1:3) * hubble + 50.0) * posunits
      
      ! Set vel in units of km/s
      vel(1:3,i) = temp(5:7)

      ! Set mass in units of 10^{10} Msun/h (Gadget convention)
      mass(i) = 10**(temp(4)-10.0) * hubble
     
      ! Set id of particle
      id(i) = i
    end do
    close(lin)

    ! Verbose some useful quantities
    write(*,*) 'Min x y z / Max x y z:'
    write(*,*) minval(pos(1,:)), minval(pos(2,:)), minval(pos(3,:))
    write(*,*) maxval(pos(1,:)), maxval(pos(2,:)), maxval(pos(3,:))
    write(*,*) 'Mass_min = ', minval(mass(:)), ' 10^10 Msun/h'
    write(*,*) 'Mass_max = ', maxval(mass(:)), ' 10^10 Msun/h'

    ! Write gadget-file
    call write_gadget(fileout,0,1,nparttotal,pos,vel,mass,id) 
 
    ! Clean up memory
    deallocate(pos)
    deallocate(vel)
    deallocate(mass)
    deallocate(id)

  end subroutine gadget_from_ascii

end module write_gadget_files

program convert
  use write_gadget_files
  implicit none

  call gadget_from_ascii()

end program convert

