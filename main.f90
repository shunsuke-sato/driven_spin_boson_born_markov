module global_variables
  implicit none
! math parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! quantum system
  complex(8) :: zrho_dm(2,2), zU_prop(2,2)
  complex(8),allocatable :: zrho_dm_memory(:,:,:),zu_prop_memory(:,:,:)
  
! time propagation
  integer :: nt, nt_memory
  real(8) :: dt, Tprop, T_memory_cut

! field parameters
  real(8) :: E0, omega0

! bath parameter  
  complex(8),allocatable :: zcorr_bath(:)
  real(8) :: omega_c, eta, beta_temp
  
  
end module global_variables
!--------------------------------------------------------------------------------------
program main
  implicit none
  use global_variables

  call input  
  call initialize

  call propagation
  
end program main
!--------------------------------------------------------------------------------------
subroutine input
  implicit none
  use global_variables


  Tprop = 200d0
  dt = 0.1d0


! laser  
  E0 = 0.1d0
  omega0 = 1d0

! bath  
  omega_c = 1d0
  eta = 1d0
  beta_temp = 1d0
  T_memory_cut  = 10d0

  
  nt = aint(Tprop/dt) + 1
  nt_memory = aint(T_memory_cut/dt) + 1



  
  
end subroutine input
!--------------------------------------------------------------------------------------
subroutine initialization
  implicit none
  use global_variables
  real(8) :: tt
  integer :: it

  allocate(zrho_dm_memory(2,2,0:nt_memory),zu_prop_memory(2,2,0:nt_memory))
  allocate(zcorr_bath(0:nt_memory))

  zrho_dm = 0d0
  zrho_dm(2,2) = 1d0

  zu_prop = 0d0
  zu_prp(1,1) = 1d0; zu_prp(2,2) = 1d0

  do it = 0, nt_memory
     tt = dt*it
! high-temperature correlation function     
     zcorr_bath(it) = 2d0*omega_c*eta/beta_temp/(1d0+(omega_c*tt)**2)
  end do
  zcorr_bath(0) = 0.5d0*zcorr_bath(0)

  zrho_dm_memory = 0d0
  zrho_dm_memory(:,:,0) = 0.5d0*zrho_dm(:,:)
  
  zu_prop_memory = 0d0
  zu_prop_memory(:,:,0) = 0.5d0*zu_prop(:,:)
  
  
  
end subroutine initialization
!--------------------------------------------------------------------------------------
subroutine propagation
  implicit none
  use global_variables
  integer :: it

  do it = 0, nt

     call dt_evolve(it)
     
  end do
  
end subroutine propagation
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
