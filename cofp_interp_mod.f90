!!!
! G.K.H. Lee - Sep 2020
! A self contained module to interpolate JANAF tables to find T and composition
! dependent heat capaicities c_p(T) for use in RT as well as future c_p(T) dynamics
! An example program is given below the module, comment this out when compiling
! with other code.
!! TODO: read mol_w directly from a file to avoid passing it.
!!!

module cofp_interp_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128


  real(kind=dp), parameter :: R = 8.31446261815324_dp


  real(kind=dp) :: cp_bar, R_bar, kappa_prime

  ! Temperature grid has 61 points from 100-6000 K in JANAF tables
  ! arrays for T and cp grid
  integer :: n_lines = 61
  real(kind=dp), allocatable, dimension(:) :: T
  real(kind=dp), allocatable, dimension(:,:) :: cp


  logical :: first_call = .True., restart = .False.

  public :: cofp_interp, cp_bar, R_bar, kappa_prime
  private :: first_call, init_cofp_interp, locate, linear_log_interp


contains

  subroutine cofp_interp(n_mol, m_list, mol_w, VMR, T_in, mmw, restart)
    implicit none

    !! Input arrays:
    ! - n_mol = number of species to compute
    ! - m_list = string array of species to compute
    ! - mol_w = array of molecules molecular weights [g mol-1] (Same order as m_list)
    ! - VMR = Volume mixing ratio of each species (Same order as m_list)
    ! - T_in = input temperature [K]
    ! - mmw - mean molecular weight (g mol-1)
    ! - restart = optional logical to restart the initialisation process (useful for testing)
    integer, intent(in) :: n_mol
    character(len=4), dimension(n_mol), intent(in) :: m_list
    real(kind=dp), intent(in) :: T_in, mmw
    real(kind=dp), dimension(n_mol), intent(in) :: mol_w, VMR
    logical, optional :: restart

    !! Work arrays
    integer :: n, iT, iT1
    real(kind=dp) :: mmr, cp_val, x0, x1, y0, y1

    !! On first call, initalise arrays and read JANAF table data for each species
    if ((first_call .eqv. .True.) .or. (restart .eqv. .True.)) then
      call init_cofp_interp(n_mol, m_list, restart)
      first_call = .False.
      restart = .False.
    end if

    !! Locate index for temperature interpolation - will need to add checks if
    ! interested in T_in < 100 K or T_in > 6000 K !!!
    call locate(T, T_in, iT)
    iT1 = iT + 1

    !! Interpolation values
    x0 = T(iT)
    x1 = T(iT1)

    !! Main loop calculations
    R_bar = 0.0_dp
    cp_bar = 0.0_dp
    do n = 1, n_mol
      ! Mass mixing ratio = VMR * molecular weight / mean molecular weight
      mmr = VMR(n) * mol_w(n)/mmw
      ! Contribution to R_bar = mmr * specific gas constant
      R_bar = R_bar + mmr * R / mol_w(n)

      ! y values for interpolation
      y0 = cp(n,iT)
      y1 = cp(n,iT1)
      call linear_log_interp(T_in,x0,x1,y0,y1,cp_val)

      ! Contribution to cp_bar = mmr * cp_val / molecular weight
      cp_bar = cp_bar + mmr * cp_val/mol_w(n)

    end do

    ! Convert R_bar [J g-1 K-1] to SI units [J kg-1 K-1]
    R_bar = R_bar * 1000.0_dp

    ! Convert cp_bar [J g-1 K-1] to SI units [J kg-1 K-1]
    cp_bar = cp_bar * 1000.0_dp

    ! kappa_prime evaluation
    kappa_prime = R_bar / cp_bar

  end subroutine cofp_interp

  subroutine init_cofp_interp(n_mol, m_list, restart)
    implicit none

    integer, intent(in) :: n_mol
    character(len=4), dimension(n_mol), intent(in) :: m_list
    logical, optional, intent(in) :: restart

    integer :: n, u, l

    ! If restart is .True. then deallocate arrays
    if (restart .eqv. .True.) then
      deallocate(T, cp)
    end if

    ! Allocate T and cp array
    allocate(T(n_lines), cp(n_mol,n_lines))

    ! Read in JANAF data for each species in order
    do n = 1, n_mol

      !print*, trim(m_list(n))

      open(newunit=u, file='data/'//trim(m_list(n))//'_JANAF.txt',action='read')
      read(u,*); read(u,*); read(u,*) ! Read 3 line header

      do l = 1, n_lines
        read(u,*) T(l), cp(n,l)
        !print*, l, T(l), cp(n,l)
      end do

      close(u)

    end do

  end subroutine init_cofp_interp

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(kind=dp), intent(in) :: xval, y1, y2, x1, x2
    real(kind=dp) :: lxval, ly1, ly2, lx1, lx2
    real(kind=dp), intent(out) :: yval
    real(kind=dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    if (lx1 >= lx2) then
      print*, 'Error in linear_log_interp: lx1 >= lx2 - STOPPING', lx1, lx2
      stop
    end if

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp


  pure subroutine locate(arr, var, idx)
    implicit none

    integer, intent(out) :: idx
    real(kind=dp), dimension(:), intent(in) :: arr
    real(kind=dp),intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section/binary search (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = size(arr)+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if (var > arr(jm)) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate


end module cofp_interp_mod


!!!
! Test program - comment out section below when coupling
!!!


program test_cofp_interp
  use cofp_interp_mod, only : cofp_interp, cp_bar, R_bar, kappa_prime
  implicit none

  integer, parameter :: n_mol = 5
  character(len=4), dimension(n_mol) :: m_list
  double precision :: T_in, mmw
  double precision, dimension(n_mol) :: mol_w, VMR
  logical :: restart

  ! Test for Earth air - gets very close to the wikipedia quoted values for cp air and R air
  m_list = (/'N2  ','O2  ','Ar  ', 'CO2 ', 'H2O '/)
  mol_w = (/28.01340d0,31.99880d0,39.9480d0,44.0095d0,18.01528d0/)
  VMR = (/0.78d0,0.21d0,0.0093d0,365d-6,1d-4/)
  T_in = 298.15d0
  mmw = 28.9645d0
  restart = .False.

  call cofp_interp(n_mol, m_list, mol_w, VMR, T_in, mmw, restart)

  print*, cp_bar, R_bar, kappa_prime, 2d0/7d0

end program test_cofp_interp
