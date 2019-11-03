module global_variables
  implicit none
! math parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! quantum system
  complex(8) :: zrho_dm(2,2), zU_prop(2,2)
  complex(8),allocatable :: zrho_dm_memory(:,:,:),zu_prop_memory(:,:,:)
  complex(8),allocatable :: zSz_memory(:,:,:)
  
! time propagation
  integer :: nt, nt_memory
  real(8) :: dt, Tprop, T_memory_cut

! field parameters
  real(8) :: E0, omega0

! bath parameter  
  complex(8),allocatable :: zcorr_bath(:)
  real(8) :: omega_c, eta, beta_temp

! Pauli matrix
  complex(8) :: zSx(2,2),zSz(2,2)
  
  
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
  allocate(zAt_memory(2,2,0:nt_memory))
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

  zSx = 0d0
  zSx(0,1) = 1d0; zSx(1,0) = 1d0 
  zSz = 0d0
  zSz(1,1) = 1d0; zSz(2,2) = -1d0

  zAt_memory(:,:,0) = zSz
  
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
subroutine dt_evolve(it)
  implicit none
  use global_variables
  integer,intent(in) :: it
  integer :: it_t
  complex(8) :: z_drho_dt(2,2), z_drho_dt_pred(2,2)
  complex(8) :: zrho_dm_old(2,2)
  complex(8) :: zA_t(2,2)

  zrho_dm_old = zrho_dm
  
  if(it /= 0)then
     z_drho_dt = 0d0
  else
     zA_t = matmul(transpose(conjg(zu_prop)), matmul(zSz, zu_prop))
     z_drho_dt = 0d0
     do it_t = 0, min(it, nt_memory)

        z_drho_dt = z_drho_dt + zcorr_bath(it_t)*(&
             matmul(zA_t,matmul(zAt_memory(:,:,it_t),zrho_dm_memory(:,:,it_t))) &
            -matmul(matmul(zAt_memory(:,:,it_t),zrho_dm_memory(:,:,it_t)),zA_t)
     end do
     zA_t = -z_drho_dt*dt
     z_drho_dt = zA_t + transpose(conjg(zA_t))
  end if

! predictor       
  z_drho_dt_pred = z_drho_dt
     
  zrho_dm = zrho_dm + dt*z_drho_dt
  call dt_evolve_zu_prop(it)
  zA_t = matmul(transpose(conjg(zu_prop)), matmul(zSz, zu_prop))
  
  do it_t = 0, nt_memory-1
     zrho_dm_memory(:,:,nt_memory-it_t) = zrho_dm_memory(:,:,nt_memory-it_t-1)
     zu_prop_memory(:,:,nt_memory-it_t) = zu_prop_memory(:,:,nt_memory-it_t-1)
     zAt_memory(:,:,nt_memory-it_t) = zAt_memory(:,:,nt_memory-it_t-1)
  end do
  zrho_dm_memory(:,:,0) = zrho_dm(:,:)
  zu_prop_memory(:,:,0) = zu_prop(:,:)
  zAt_memory(:,:,0) = zA_t(:,:)

! corrector     
  z_drho_dt = 0d0
  do it_t = 0, min(it+1, nt_memory)
     
     z_drho_dt = z_drho_dt + zcorr_bath(it_t)*(&
          matmul(zA_t,matmul(zAt_memory(:,:,it_t),zrho_dm_memory(:,:,it_t))) &
          -matmul(matmul(zAt_memory(:,:,it_t),zrho_dm_memory(:,:,it_t)),zA_t)
  end do
  zA_t = -z_drho_dt*dt
  z_drho_dt = zA_t + transpose(conjg(zA_t))

  zrho_dm = zrho_dm_old + dt*0.5d0*(z_drho_dt + z_drho_dt_pred)
  zrho_dm_memory(:,:,0) = zrho_dm(:,:)  

     
end subroutine dt_evolve
!--------------------------------------------------------------------------------------
subroutine dt_evolve_zu_prop(it)
  implicit none
  use global_variables
  integer,intent(in) :: it
  complex(8) :: zpropagator(2,2),zpropagator_t(2,2)
  real(8) :: Et, Ham(2,2), tt
  real(8) :: eig_vec(2,2), eig_val(2)
  real(8) :: ss

  tt = dt*it+0.5d0*dt
  Et = E0*sin(omega0*tt)

  Ham(1,1) =  0.5d0
  Ham(2,2) = -0.5d0
  Ham(1,2) = Et
  Ham(2,1) = Et

  eig_val(1) = 0.5d0*(ham(1,1)+ham(2,2)+sqrt((ham(2,2)-ham(1,1))**2+4d0*ham(1,2)**2))
  eig_val(2) = 0.5d0*(ham(1,1)+ham(2,2)-sqrt((ham(2,2)-ham(1,1))**2+4d0*ham(1,2)**2))
  eig_vec(1,1) = 1d0; eig_vec(2,1) = ham(1,2)/(eig_val(1)-ham(2,2))
  eig_vec(1,2) = ham(1,2)/(eig_val(2)-ham(1,1)); eig_vec(2,2) = 1d0

  ss = eig_vec(1,1)**2 + eig_vec(2,1)**2
  eig_vec(:,1) = eig_vec(:,1)/sqrt(ss)

  ss = eig_vec(1,2)**2 + eig_vec(2,2)**2
  eig_vec(:,2) = eig_vec(:,2)/sqrt(ss)

  zpropagator = 0d0
  zpropagator(1,1) = exp(-zi*dt*eig_val(1))
  zpropagator(2,2) = exp(-zi*dt*eig_val(2))

  zpropagator_t = matmul(eig_val, matmul(zpropagator,transpose(eig_vec)))

  zu_prop = matmul(zpropagator_t,zu_prop)
  
  
  

end subroutine dt_evolve_zu_prop
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
