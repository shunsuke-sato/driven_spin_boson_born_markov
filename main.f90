module global_variables
  implicit none
! math parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! quantum system
  complex(8) :: zrho_dm(2,2), zU_prop(2,2), zrho_dm_s(2,2)
  complex(8),allocatable :: zrho_dm_memory(:,:,:),zu_prop_memory(:,:,:)
  complex(8),allocatable :: zAt_memory(:,:,:),zAt_zrho_memory(:,:,:)
  
! time propagation
  integer :: nt
  real(8) :: dt, Tprop, T_memory_cut

! field parameters
  real(8) :: E0, omega0

! bath parameter  
  complex(8),allocatable :: zcorr_bath(:)
  real(8) :: omega_c, eta, beta_temp
  real(8) :: rx, rz
  complex(8) :: zGamma_bath_zero, zGamma_bath_p_delta, zGamma_bath_m_delta
  real(8) :: gamma_bath_zero, gamma_bath_p_delta, gamma_bath_m_delta
  real(8) :: s_bath_zero, s_bath_p_delta, s_bath_m_delta
  

! Pauli matrix
  complex(8) :: zSx(2,2),zSy(2,2),zSz(2,2)
  

! propagation scheme
  integer :: n_propagation_scheme
  integer,parameter :: N_propagation_Born = 0
  integer,parameter :: N_propagation_Redfield = 1
  integer,parameter :: N_propagation_Lindblad = 2

  integer,parameter :: N_propagation_fidelity_analysis = -1

! lindblad eq.
  real(8) :: H_LS(2,2), lamb_shift(2), rho_eq(2)
  real(8) :: T1, T2

! Floquet analysis
  logical,parameter :: if_floquet_analysis = .true.
  integer,parameter :: nf_cut = 30
  integer :: nt_floquet_cycle
  complex(8),allocatable :: zstates_floquet(:,:)
  real(8),allocatable :: eps_floquet(:)
  
  
end module global_variables
!--------------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input  

  select case(n_propagation_scheme)
    case(N_propagation_Born,N_propagation_Redfield)
      call initialization
      call propagation
    case(N_propagation_Lindblad)
      call initialization_Lindblad
      call propagation_Lindblad
    case(N_propagation_fidelity_analysis)
      call fidelity_analysis
    case default
      stop 'Error: Invalid propagation scheme'
    end select
  
end program main
!--------------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none
  real(8) :: ancycle

! scheme
!  n_propagation_scheme = N_propagation_Born
!  n_propagation_scheme = N_propagation_Lindblad
  n_propagation_scheme = N_propagation_fidelity_analysis

! laser  
  E0 = 0.50d0
  omega0 = 1.0d0

! propagation
  Tprop = 400d0
  dt = 0.01d0

! bath  
  omega_c = 0.5d0
  eta = 0.1d0 !1d0 ! debug
  beta_temp = 1d0
  T_memory_cut  = 10d0
  rx = 1d0
  rz = 0d0

! optimizing parameters
  write(*,"(A,2x,e16.6e3)")"Tprop (input) =", Tprop

  ancycle = aint(max(Tprop, 0d0)/(2d0*pi/omega0))
  Tprop = (1+ancycle)*2d0*pi/omega0

  write(*,"(A,2x,e26.16e3)")"Tprop (mod.)  =", Tprop

  write(*,"(A,2x,e26.16e3)")"dt (input) =", dt
  nt = nint(Tprop/dt) + 1
  dt = Tprop/nt

  nt_floquet_cycle = nint((2d0*pi/omega0)/dt)
  write(*,*)'aint',aint(Tprop/dt)
  write(*,"(A,2x,e26.16e3)")"dt (mod.) =", dt
  write(*,"(A,2x,I9)")"nt (mod.) =", nt

  zSx = 0d0
  zSx(1,2) = 1d0; zSx(2,1) = 1d0 
  zSy = 0d0
  zSy(1,2) = -zi; zSy(2,1) = zi
  zSz = 0d0
  zSz(1,1) = 1d0; zSz(2,2) = -1d0
  
end subroutine input
!--------------------------------------------------------------------------------------
subroutine initialization
  use global_variables
  implicit none

  allocate(zrho_dm_memory(2,2,0:nt+1),zu_prop_memory(2,2,0:nt+1))
  allocate(zAt_memory(2,2,0:nt+1),zAt_zrho_memory(2,2,0:nt+1))

  zrho_dm = 0d0
  zrho_dm(2,2) = 1d0

  zu_prop = 0d0
  zu_prop(1,1) = 1d0; zu_prop(2,2) = 1d0

  call calc_bath_correlation



  zrho_dm_memory = 0d0
  zrho_dm_memory(:,:,0) = zrho_dm(:,:)
  
  zu_prop_memory = 0d0
  zu_prop_memory(:,:,0) = zu_prop(:,:)

  if(if_floquet_analysis) call calc_Floquet_states

  
end subroutine initialization
!--------------------------------------------------------------------------------------
subroutine initialization_Lindblad
  use global_variables
  implicit none
  real(8) :: tt
  integer :: it


  zrho_dm = 0d0
  zrho_dm(2,2) = 1d0

  zrho_dm_s = zrho_dm

  call calc_bath_correlation
  call FT_bath_correlation


  lamb_shift(1) = s_bath_zero*rz**2 + s_bath_p_delta*rx**2
  lamb_shift(2) = s_bath_zero*rz**2 + s_bath_m_delta*rx**2


  rho_eq(1) = gamma_bath_m_delta/(gamma_bath_p_delta+gamma_bath_m_delta)
  rho_eq(2) = gamma_bath_p_delta/(gamma_bath_p_delta+gamma_bath_m_delta)

  T1 = rx**2*(gamma_bath_p_delta+gamma_bath_m_delta)
  T1 = 1d0/T1

  T2 = 0.5d0*rx**2*(gamma_bath_p_delta+gamma_bath_m_delta) + rz**2*gamma_bath_zero
  T2 = 1d0/T2

  write(*,"(A,2x,999e26.16e3)")"lamb_shift=",lamb_shift
  write(*,"(A,2x,999e26.16e3)")"T1=",T1
  write(*,"(A,2x,999e26.16e3)")"T2=",T2
  if(if_floquet_analysis) call calc_Floquet_states
  
end subroutine initialization_Lindblad
!--------------------------------------------------------------------------------------
subroutine calc_bath_correlation
  use global_variables
  implicit none
  real(8) :: tt
  integer :: it

  allocate(zcorr_bath(0:nt+1))
  do it = 0, nt+1
     tt = dt*it
! high-temperature correlation function     
!     zcorr_bath(it) = 2d0*omega_c*eta/beta_temp/(1d0+(omega_c*tt)**2)
! zero-temperature correlation function     
     zcorr_bath(it) = eta*omega_c**2/(1d0+zi*omega_c*tt)**2
  end do

end subroutine calc_bath_correlation
!--------------------------------------------------------------------------------------
subroutine FT_bath_correlation
  use global_variables
  implicit none
  real(8) :: tt
  integer :: it
  

  zGamma_bath_zero = 0d0
  zGamma_bath_p_delta = 0d0
  zGamma_bath_m_delta = 0d0
  do it = 0, nt+1
    tt = dt*it

    zGamma_bath_zero = zGamma_bath_zero + zcorr_bath(it)
    zGamma_bath_p_delta = zGamma_bath_p_delta + zcorr_bath(it)*exp(zi*tt)
    zGamma_bath_m_delta = zGamma_bath_m_delta + zcorr_bath(it)*exp(-zi*tt)

  end do

  zGamma_bath_zero = zGamma_bath_zero*dt
  zGamma_bath_p_delta = zGamma_bath_p_delta*dt
  zGamma_bath_m_delta = zGamma_bath_m_delta*dt

  gamma_bath_zero = 2d0*real(zGamma_bath_zero)
  gamma_bath_p_delta= 2d0*real(zGamma_bath_p_delta)
  gamma_bath_m_delta= 2d0*real(zGamma_bath_m_delta)

  s_bath_zero = aimag(zGamma_bath_zero)
  s_bath_p_delta = aimag(zGamma_bath_p_delta)
  s_bath_m_delta = aimag(zGamma_bath_m_delta)


end subroutine FT_bath_correlation
!--------------------------------------------------------------------------------------
subroutine propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: S_F_fidelity, S_F_fidelity_ave
  character(256) :: cmethod

  if(n_propagation_scheme == N_propagation_Born)then
    cmethod='born'
  else if(n_propagation_scheme == N_propagation_Redfield)then
    cmethod='redfield'
  end if


  if(if_floquet_analysis)open(31,file="floquet_fidelity_"//trim(cmethod)//".out")
  S_F_fidelity_ave = 0d0
  call pre_propagation

  open(20,file='pop_t_'//trim(cmethod)//'.out')
  do it = 0, nt

    zrho_dm_s = matmul(zu_prop_memory(:,:,it), &
      matmul(zrho_dm, conjg(transpose(zu_prop_memory(:,:,it)))))
    write(20,"(999e26.16e3)")dt*it,real(zrho_dm_s(1,1)),real(zrho_dm_s(2,2)),zrho_dm_s(1,2)

    if(if_floquet_analysis .and. it >= nt-nt_floquet_cycle+1)then
      call calc_instantaneous_floquet_fidelity(zrho_dm_s, S_F_fidelity, it*dt)
      write(31,"(999e26.16e3)")dt*it,S_F_fidelity
      S_F_fidelity_ave = S_F_fidelity_ave + S_F_fidelity
    end if

    if(n_propagation_scheme == N_propagation_Born)then
      call dt_evolve(it)
    else if(n_propagation_scheme == N_propagation_Redfield)then
      call dt_evolve_Redfield(it)
    else
      stop 'Error in propagation'
    end if
     
  end do
  close(20)


  if(if_floquet_analysis)then
    write(*,"(A,2x,e16.6e3)")'Floquet fidelity (cycle averaged)=',S_F_fidelity_ave/nt_floquet_cycle

    close(31)
  end if
end subroutine propagation
!--------------------------------------------------------------------------------------
subroutine propagation_lindblad
  use global_variables
  implicit none
  integer :: it
  real(8) :: S_F_fidelity, S_F_fidelity_ave


  if(if_floquet_analysis)open(31,file="floquet_fidelity_ldndblad.out")
  S_F_fidelity_ave = 0d0
  open(20,file='pop_t_lindblad.out')
  do it = 0, nt

    write(20,"(999e26.16e3)")dt*it,real(zrho_dm_s(1,1)),real(zrho_dm_s(2,2)),zrho_dm_s(1,2)

    if(if_floquet_analysis .and. it >= nt-nt_floquet_cycle+1)then
      call calc_instantaneous_floquet_fidelity(zrho_dm_s, S_F_fidelity, it*dt)
      write(31,"(999e26.16e3)")dt*it,S_F_fidelity
      S_F_fidelity_ave = S_F_fidelity_ave + S_F_fidelity
    end if
    call dt_evolve_lindblad(it)
     
  end do
  close(20)

  if(if_floquet_analysis)then
    write(*,"(A,2x,e16.6e3)")'Floquet fidelity (cycle averaged)=',S_F_fidelity_ave/nt_floquet_cycle

    close(31)
  end if
end subroutine propagation_lindblad
!--------------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: it_t
  complex(8) :: z_drho_dt(2,2), z_drho_dt_pred(2,2)
  complex(8) :: zrho_dm_old(2,2)
  complex(8) :: zA_t(2,2)

  zrho_dm_old = zrho_dm


! predictor
  if(it/=0)then
    z_drho_dt(:,:) = -0.5d0*zcorr_bath(0) &
      *(matmul(zAt_memory(:,:,it),zAt_zrho_memory(:,:,it)) &
      -matmul(zAt_zrho_memory(:,:,it),zAt_memory(:,:,it)))

    do it_t = 1, it-1
      z_drho_dt(:,:) = z_drho_dt(:,:) &
        -zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it),zAt_zrho_memory(:,:,it-it_t)) &
        -matmul(zAt_zrho_memory(:,:,it-it_t),zAt_memory(:,:,it)))
    end do

    it_t = it
    z_drho_dt(:,:) = z_drho_dt(:,:) &
      -0.5d0*zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it),zAt_zrho_memory(:,:,it-it_t)) &
      -matmul(zAt_zrho_memory(:,:,it-it_t),zAt_memory(:,:,it)))
  else
    z_drho_dt = 0d0
  end if

  z_drho_dt(:,:) = z_drho_dt(:,:) +  transpose(conjg(z_drho_dt(:,:)))
  z_drho_dt = z_drho_dt*dt

  z_drho_dt_pred = z_drho_dt

  zrho_dm = zrho_dm + dt* z_drho_dt
  zAt_zrho_memory(:,:,it+1) = matmul(zAt_memory(:,:,it+1),zrho_dm(:,:))

! corrector
    z_drho_dt(:,:) = &
      -0.5d0*zcorr_bath(0) &
      *(matmul(zAt_memory(:,:,it+1),zAt_zrho_memory(:,:,it+1)) &
      -matmul(zAt_zrho_memory(:,:,it+1),zAt_memory(:,:,it+1)))

    do it_t = 1, it+1-1
      z_drho_dt(:,:) = z_drho_dt(:,:) &
        -zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it+1),zAt_zrho_memory(:,:,it+1-it_t)) &
        -matmul(zAt_zrho_memory(:,:,it+1-it_t),zAt_memory(:,:,it+1)))
    end do

    it_t = it+1
    z_drho_dt(:,:) = z_drho_dt(:,:) &
      -0.5d0*zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it+1),zAt_zrho_memory(:,:,it+1-it_t)) &
      -matmul(zAt_zrho_memory(:,:,it+1-it_t),zAt_memory(:,:,it+1)))

    z_drho_dt(:,:) = z_drho_dt(:,:) +  transpose(conjg(z_drho_dt(:,:)))
    z_drho_dt = z_drho_dt*dt

    zrho_dm = zrho_dm_old + 0.5d0*dt*(z_drho_dt + z_drho_dt_pred)
    zrho_dm_memory(:,:,it+1) = zrho_dm(:,:)
    zAt_zrho_memory(:,:,it+1) = matmul(zAt_memory(:,:,it+1),zrho_dm(:,:))


     
end subroutine dt_evolve
!--------------------------------------------------------------------------------------
subroutine dt_evolve_Redfield(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: it_t
  complex(8) :: z_drho_dt(2,2), z_drho_dt_pred(2,2)
  complex(8) :: zrho_dm_old(2,2)
  complex(8) :: zA_t(2,2), zAt_zrhot(2,2)

  zrho_dm_old = zrho_dm


! predictor
  if(it/=0)then

    it_t = 0
    zAt_zrhot(:,:) = matmul(zAt_memory(:,:, it-it_t), zrho_dm)
    z_drho_dt(:,:) = -0.5d0*zcorr_bath(0) &
      *(matmul(zAt_memory(:,:,it),zAt_zrhot(:,:)) &
      -matmul(zAt_zrhot(:,:),zAt_memory(:,:,it)))

    do it_t = 1, it-1
      zAt_zrhot(:,:) = matmul(zAt_memory(:,:, it-it_t), zrho_dm)
      z_drho_dt(:,:) = z_drho_dt(:,:) &
        -zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it),zAt_zrhot(:,:)) &
        -matmul(zAt_zrhot(:,:),zAt_memory(:,:,it)))
    end do

    it_t = it
    zAt_zrhot(:,:) = matmul(zAt_memory(:,:, it-it_t), zrho_dm)
    z_drho_dt(:,:) = z_drho_dt(:,:) &
      -0.5d0*zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it),zAt_zrhot(:,:)) &
      -matmul(zAt_zrhot(:,:),zAt_memory(:,:,it)))
  else
    z_drho_dt = 0d0
  end if

  z_drho_dt(:,:) = z_drho_dt(:,:) +  transpose(conjg(z_drho_dt(:,:)))
  z_drho_dt = z_drho_dt*dt

  z_drho_dt_pred = z_drho_dt

  zrho_dm = zrho_dm + dt* z_drho_dt

! corrector
    it_t = 0
    zAt_zrhot(:,:) = matmul(zAt_memory(:,:, it+1-it_t), zrho_dm)

    z_drho_dt(:,:) = &
      -0.5d0*zcorr_bath(0) &
      *(matmul(zAt_memory(:,:,it+1),zAt_zrhot(:,:)) &
      -matmul(zAt_zrhot(:,:),zAt_memory(:,:,it+1)))

    do it_t = 1, it+1-1
      zAt_zrhot(:,:) = matmul(zAt_memory(:,:, it+1-it_t), zrho_dm)
      z_drho_dt(:,:) = z_drho_dt(:,:) &
        -zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it+1),zAt_zrhot(:,:)) &
        -matmul(zAt_zrhot(:,:),zAt_memory(:,:,it+1)))
    end do

    it_t = it+1
    zAt_zrhot(:,:) = matmul(zAt_memory(:,:, it+1-it_t), zrho_dm)
    z_drho_dt(:,:) = z_drho_dt(:,:) &
      -0.5d0*zcorr_bath(it_t)*(matmul(zAt_memory(:,:,it+1),zAt_zrhot(:,:)) &
      -matmul(zAt_zrhot(:,:),zAt_memory(:,:,it+1)))

    z_drho_dt(:,:) = z_drho_dt(:,:) +  transpose(conjg(z_drho_dt(:,:)))
    z_drho_dt = z_drho_dt*dt

    zrho_dm = zrho_dm_old + 0.5d0*dt*(z_drho_dt + z_drho_dt_pred)

     
end subroutine dt_evolve_Redfield
!--------------------------------------------------------------------------------------
subroutine dt_evolve_lindblad(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  complex(8) :: k1(2,2), k2(2,2), k3(2,2), k4(2,2), zrho_tmp(2,2)
  real(8) :: Et, Ham(2,2), tt
  real(8) :: ss


! k1
  zrho_tmp = zrho_dm_s
  tt = dt*it

  call calc_ham_matrix(ham,tt)
  Ham(1,1) = Ham(1,1) + lamb_shift(1)
  Ham(2,2) = Ham(2,2) + lamb_shift(2)

  k1 = -zi*(matmul(ham,zrho_tmp)-matmul(zrho_tmp,ham))
  k1(1,1) = k1(1,1) -(1d0/T1)*(zrho_tmp(1,1)-rho_eq(1))
  k1(2,1) = k1(2,1) -(1d0/T2)*zrho_tmp(2,1)
  k1(1,2) = k1(1,2) -(1d0/T2)*zrho_tmp(1,2)
  k1(2,2) = k1(2,2) -(1d0/T1)*(zrho_tmp(2,2)-rho_eq(2))

! k2
  zrho_tmp = zrho_dm_s + 0.5d0*dt*k1
  tt = dt*it+0.5d0*dt

  call calc_ham_matrix(ham,tt)
  Ham(1,1) = Ham(1,1) + lamb_shift(1)
  Ham(2,2) = Ham(2,2) + lamb_shift(2)

  k2 = -zi*(matmul(ham,zrho_tmp)-matmul(zrho_tmp,ham))
  k2(1,1) = k2(1,1) -(1d0/T1)*(zrho_tmp(1,1)-rho_eq(1))
  k2(2,1) = k2(2,1) -(1d0/T2)*zrho_tmp(2,1)
  k2(1,2) = k2(1,2) -(1d0/T2)*zrho_tmp(1,2)
  k2(2,2) = k2(2,2) -(1d0/T1)*(zrho_tmp(2,2)-rho_eq(2))

! k3
  zrho_tmp = zrho_dm_s + 0.5d0*dt*k2
  tt = dt*it+0.5d0*dt

  call calc_ham_matrix(ham,tt)
  Ham(1,1) = Ham(1,1) + lamb_shift(1)
  Ham(2,2) = Ham(2,2) + lamb_shift(2)

  k3 = -zi*(matmul(ham,zrho_tmp)-matmul(zrho_tmp,ham))
  k3(1,1) = k3(1,1) -(1d0/T1)*(zrho_tmp(1,1)-rho_eq(1))
  k3(2,1) = k3(2,1) -(1d0/T2)*zrho_tmp(2,1)
  k3(1,2) = k3(1,2) -(1d0/T2)*zrho_tmp(1,2)
  k3(2,2) = k3(2,2) -(1d0/T1)*(zrho_tmp(2,2)-rho_eq(2))

! k4
  zrho_tmp = zrho_dm_s + dt*k3
  tt = dt*it+dt

  call calc_ham_matrix(ham,tt)
  Ham(1,1) = Ham(1,1) + lamb_shift(1)
  Ham(2,2) = Ham(2,2) + lamb_shift(2)

  k4 = -zi*(matmul(ham,zrho_tmp)-matmul(zrho_tmp,ham))
  k4(1,1) = k4(1,1) -(1d0/T1)*(zrho_tmp(1,1)-rho_eq(1))
  k4(2,1) = k4(2,1) -(1d0/T2)*zrho_tmp(2,1)
  k4(1,2) = k4(1,2) -(1d0/T2)*zrho_tmp(1,2)
  k4(2,2) = k4(2,2) -(1d0/T1)*(zrho_tmp(2,2)-rho_eq(2))


! sum
  
  zrho_dm_s = zrho_dm_s + (dt/6d0)*(k1+2d0*k2+2d0*k3+k4)

  

end subroutine dt_evolve_lindblad
!--------------------------------------------------------------------------------------
subroutine dt_evolve_zu_prop(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  complex(8) :: zpropagator(2,2),zpropagator_t(2,2)
  real(8) :: Et, Ham(2,2), tt
  real(8) :: eig_vec(2,2), eig_val(2)
  real(8) :: ss

  tt = dt*it+0.5d0*dt
  call calc_ham_matrix(ham,tt)

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

  zpropagator_t = matmul(eig_vec, matmul(zpropagator,transpose(eig_vec)))

  zu_prop = matmul(zpropagator_t,zu_prop)
  
  
  

end subroutine dt_evolve_zu_prop
!--------------------------------------------------------------------------------------
subroutine calc_ham_matrix(ham_mat, tt_in)
  use global_variables
  implicit none
  real(8),intent(out) :: ham_mat(2,2)
  real(8),intent(in) :: tt_in
  real(8) :: Et

  Et = E0*sin(omega0*tt_in)

  ham_mat(1,1) =  0.5d0
  ham_mat(2,2) = -0.5d0
  ham_mat(1,2) = Et
  ham_mat(2,1) = Et

end subroutine calc_ham_matrix
!--------------------------------------------------------------------------------------
subroutine pre_propagation
  use global_variables
  implicit none
  integer :: it

! propagator
  do it = 0, nt

    call dt_evolve_zu_prop(it)
! here zu_prop is the subsystem propagator from 0 to (it+1)*dt
    zu_prop_memory(:,:,it+1) = zu_prop(:,:)
    
  end do


! interaction matrix
  zAt_memory = 0d0
  do it = 0, nt+1
    zAt_memory(:,:,it) = matmul(conjg(transpose(zu_prop_memory(:,:,it))), &
      matmul(rx*zSx + rz*zSz, zu_prop_memory(:,:,it)))
  end do

  zAt_zrho_memory(:,:,0) = matmul(zAt_memory(:,:,0),zrho_dm_memory(:,:,0))

end subroutine pre_propagation
!--------------------------------------------------------------------------------------
! Here, H_ext(t) = E_0*sin(omega_0*t)*S_x is assumed
!--------------------------------------------------------------------------------------
subroutine calc_Floquet_states
  use global_variables
  implicit none
  integer :: ndim_s, ndim_e
  integer :: icut, i1, i2
  complex(8),allocatable :: zham_f(:,:)
  real(8),allocatable :: eps_f(:)
! lapack
  integer :: ndim, lwork, infor
  complex(8),allocatable :: work(:)
  real(8),allocatable :: rwork(:)

  ndim_s = 1 -2*nf_cut
  ndim_e = 2 +2*nf_cut
  allocate(zham_f(ndim_s:ndim_e,ndim_s:ndim_e))
  allocate(eps_f(ndim_s:ndim_e))
  zham_f = 0d0

  do icut = -nf_cut, nf_cut
    i1 = 1 + icut*2
    i2 = 2 + icut*2
    zham_f(i1:i2,i1:i2) = 0.5d0*zSz(1:2,1:2)
    zham_f(i1,i1) = zham_f(i1,i1) - icut*omega0
    zham_f(i2,i2) = zham_f(i2,i2) - icut*omega0

    if(icut /= nf_cut)then
      zham_f(i1+2:i2+2,i1:i2) =  0.5d0*zi*E0*zSx(1:2,1:2)
      zham_f(i1:i2,i1+2:i2+2) = -0.5d0*zi*E0*zSx(1:2,1:2)
    end if

  end do

! diagonalization
  ndim = ndim_e-ndim_s+1
  lwork = 2*(ndim+1)*ndim+64
  allocate(work(lwork), rwork(max(1,3*ndim-2)))
  call zheev('V', 'U', ndim, zham_f, ndim, eps_f, work, lwork, rwork, infor)


  allocate(zstates_floquet(ndim_s:ndim_e,ndim_s:ndim_e))
  allocate(eps_floquet(ndim_s:ndim_e))


  zstates_floquet = zham_f
  eps_floquet = eps_f
!  write(*,*)'eps_floquet(1:2)',eps_floquet(1:2)
end subroutine calc_Floquet_states
!--------------------------------------------------------------------------------------
subroutine  provide_Floquet_state_vectors_at_t(zpsi_F_out, tt_in)
  use global_variables
  implicit none
  complex(8),intent(out) :: zpsi_F_out(2,2)
  real(8),intent(in) :: tt_in
  integer :: ndim_s, ndim_e
  integer :: icut, i1, i2


  ndim_s = 1 -2*nf_cut
  ndim_e = 2 +2*nf_cut

  zpsi_F_out = 0d0
  do icut = -nf_cut, nf_cut
    i1 = 1 + icut*2
    i2 = 2 + icut*2

    zpsi_F_out(:,1) = zpsi_F_out(:,1) + exp(-zi*omega0*icut*tt_in)*zstates_floquet(i1:i2,1)
    zpsi_F_out(:,2) = zpsi_F_out(:,2) + exp(-zi*omega0*icut*tt_in)*zstates_floquet(i1:i2,2)
  end do
  zpsi_F_out(:,1) = zpsi_F_out(:,1)*exp(-zi*eps_floquet(1)*tt_in)
  zpsi_F_out(:,2) = zpsi_F_out(:,2)*exp(-zi*eps_floquet(2)*tt_in)


!  write(*,*)'norm',sum(abs(zpsi_F_out)**2), sum(conjg(zpsi_F_out(:,1))*zpsi_F_out(:,2))
  
end subroutine provide_Floquet_state_vectors_at_t
!--------------------------------------------------------------------------------------
subroutine calc_instantaneous_floquet_fidelity(zrho_in, S_F_fidelity_out, tt_in)
  use global_variables
  implicit none
  complex(8),intent(in) :: zrho_in(2,2)
  real(8),intent(in) :: tt_in
  real(8),intent(out) :: S_F_fidelity_out
  complex(8) :: zstates_nat(2,2), zpsi_F(2,2)
  real(8) :: occ_nat(2)
  real(8) :: S_F(2,2)
  integer :: i,j
  complex(8) :: zs

  S_F_fidelity_out = 0d0
  call diag_2x2(zrho_in, zstates_nat, occ_nat)
  
  call provide_Floquet_state_vectors_at_t(zpsi_F, tt_in)


  do i = 1,2
    do j = 1,2

      zs = sum(conjg(zstates_nat(:,i))*zpsi_F(:,j))
      S_F(i,j) = abs(zs)**2

    end do
  end do

!  write(*,*)S_F
  S_F_fidelity_out = abs(S_F(1,1)*S_F(2,2)-S_F(1,2)*S_F(2,1))
end subroutine calc_instantaneous_floquet_fidelity
!--------------------------------------------------------------------------------------
subroutine diag_2x2(zmat, zvec, lambda)
  implicit none
  complex(8),intent(in) :: zmat(2,2)
  complex(8),intent(out) :: zvec(2,2)
  real(8),intent(out) :: lambda(2)
  real(8) :: a, c
  complex(8) :: zb
  real(8) :: ss

  zvec = 0d0
  lambda = 0d0

  a  = zmat(1,1)
  c  = zmat(2,2)
  zb = zmat(1,2)

  lambda(1) = 0.5d0*((a+c) + sqrt((a-c)**2 + 4d0*abs(zb)**2)) 
  lambda(2) = 0.5d0*((a+c) - sqrt((a-c)**2 + 4d0*abs(zb)**2)) 


  if( abs(lambda(1) - a) > abs(lambda(1) - c)  ) then
    zvec(2,1) = 1d0
    zvec(1,1) = zb/(lambda(1)-a)

    zvec(1,2) = 1d0
    zvec(2,2) = conjg(zb)/(lambda(2)-c)
  else
    zvec(1,1) = 1d0
    zvec(2,1) = conjg(zb)/(lambda(1)-c)

    zvec(2,2) = 1d0
    zvec(1,2) = zb/(lambda(2)-a)
  end if

!  write(*,*)'Error:', sum(abs(matmul(zmat, zvec(:,1))-lambda(1)*zvec(:,1))**2) &
!      +sum(abs(matmul(zmat, zvec(:,2))-lambda(2)*zvec(:,2))**2)
  
  ss = sum(abs(zvec(:,1))**2)
  zvec(:,1) = zvec(:,1)/sqrt(ss)

  ss = sum(abs(zvec(:,2))**2)
  zvec(:,2) = zvec(:,2)/sqrt(ss)

end subroutine diag_2x2
!--------------------------------------------------------------------------------------
subroutine fidelity_analysis
  use global_variables
  implicit none
  real(8) :: a,br,bi,c, tmp
  complex(8) :: zb
  complex(8),allocatable :: zrho_dm_Born(:,:,:)
  complex(8),allocatable :: zrho_dm_Lindblad(:,:,:)
  integer :: it
  complex(8) :: zvec(2,2), zrho_dm_sqrt(2,2)
  complex(8) :: zS_fidelity_mat(2,2)
  real(8) :: occ(2), occ_sqrt(2), lambda(2)
  real(8) :: fidelity_self, fidelity_born_vs_lindblad
  

  allocate(zrho_dm_Born(2,2,0:nt))
  allocate(zrho_dm_Lindblad(2,2,0:nt))

  open(40,file='pop_t.out')
  open(41,file='pop_t_lindblad.out')
  do it = 0, nt
    read(40,*)tmp,a,c,br,bi
    zrho_dm_Born(1,1,it) = a
    zrho_dm_Born(2,1,it) = br - zi*bi
    zrho_dm_Born(1,2,it) = br + zi*bi
    zrho_dm_Born(2,2,it) = c

    read(41,*)tmp,a,c,br,bi
    zrho_dm_lindblad(1,1,it) = a
    zrho_dm_lindblad(2,1,it) = br - zi*bi
    zrho_dm_lindblad(1,2,it) = br + zi*bi
    zrho_dm_lindblad(2,2,it) = c

  end do
  close(40)
  close(41)

  open(42,file="fidelity_t.out")
! compute fidelity
  do it = 0, nt
    call diag_2x2(zrho_dm_Born(:,:,it), zvec, occ)
    occ_sqrt = sqrt(occ)
    zrho_dm_sqrt(1,1) = occ_sqrt(1)* zvec(1,1)*conjg(zvec(1,1)) &
                      + occ_sqrt(2)* zvec(1,2)*conjg(zvec(1,2))

    zrho_dm_sqrt(2,1) = occ_sqrt(1)* zvec(2,1)*conjg(zvec(1,1)) &
                      + occ_sqrt(2)* zvec(2,2)*conjg(zvec(1,2))

    zrho_dm_sqrt(1,2) = occ_sqrt(1)* zvec(1,1)*conjg(zvec(2,1)) &
                      + occ_sqrt(2)* zvec(1,2)*conjg(zvec(2,2))

    zrho_dm_sqrt(2,2) = occ_sqrt(1)* zvec(2,1)*conjg(zvec(2,1)) &
                      + occ_sqrt(2)* zvec(2,2)*conjg(zvec(2,2))


! self fidelity    

    zS_fidelity_mat = matmul(matmul(zrho_dm_sqrt, zrho_dm_born(:,:,it)),zrho_dm_sqrt)
    call diag_2x2(zS_fidelity_mat, zvec, lambda)

    fidelity_self = sum(sqrt(lambda))**2


    zS_fidelity_mat = matmul(matmul(zrho_dm_sqrt, zrho_dm_lindblad(:,:,it)),zrho_dm_sqrt)
    call diag_2x2(zS_fidelity_mat, zvec, lambda)

    fidelity_born_vs_lindblad = sum(sqrt(lambda))**2

    write(42,"(999e26.16e3)")dt*it, fidelity_self, fidelity_born_vs_lindblad

! Born vs Lindblad
  end do
  close(42)
end subroutine fidelity_analysis
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
