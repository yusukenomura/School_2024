!
!  Solving 1D antiferromagnetic Heisenberg model by RBM 
!    using stochastic reconfiguration (SR) method for optimization  
!
!  Copyright (C) 2024 Yusuke Nomura
!
!
!  We assume translational symmetry in variational parameters 
!  Bias terms a_i b_j are set to zero 
!  We only optimize interaction parameters W_ij
!
program RBM_solver 
  use mod_RBM
  implicit none 

  integer :: L         ! system size
  integer :: Nv        ! number of variational parameters 
  integer :: Nstep     ! total number of optimization steps
  integer :: Nsample   ! number of Monte Carlo samples
  integer :: Nwarmup   ! length of Warmup
  real(8) :: delta_tau ! controls the magnitude of parameter update 

  integer :: f         ! index for independent neurons 
  integer :: i         ! index for visible units
  integer :: j         ! index for hidden units
  integer :: k         ! index for variational parameters 
  integer :: jj        ! subindex for hidden units
  integer :: iw        ! index for W interaction 
  integer :: iteration ! index for optimization iterations
  integer :: iupdate   ! index for Monte Carlo update trials
  integer :: isample   ! index for Monte Carlo samples

  real(8) :: psi_x     ! psi(x)  = <x |psi>
  real(8) :: psi_xp    ! psi(x') = <x'|psi>
  real(8) :: Eloc_x    ! Eloc(x) = \sum_x'  <x|H|x'> * ( psi(x')/psi(x) )
  real(8) :: E         ! total energy 
  real(8) :: delW_max  ! maximum change in W at each optimization step

  integer , allocatable :: x (:)         ! |x>  = |sigma_1 , sigma_2 , ..., sigma_N >
  integer , allocatable :: xp(:)         ! |x'> = |sigma'_1, sigma'_2, ..., sigma'_N>
  real(8) , allocatable :: theta (:)     ! theta (j) = \sum_i  W_ij * \sigma_i
  real(8) , allocatable :: thetap(:)     ! thetap(j) = \sum_i  W_ij * \sigma'_i 
  real(8) , allocatable :: gvec(:)       ! derivative of energy with respect to variational parameters 
  real(8) , allocatable :: Smat(:,:)     ! metric  


  real(8) :: Sdiag_av
  integer :: i1, i2
  integer :: k1, k2
  integer :: iseed
  real(8) , allocatable :: Ovec(:)        
  real(8) , allocatable :: Ovec_loc_x(:)  
  real(8) , external :: grnd
  ! 
  ! set seed for random number generator  
  !
  call system_clock(count=iseed)
  call sgrnd(iseed)
  ! 
  ! read inputs
  ! 
  open(unit=1,file='RBM.input',status='old')
  read(1,*)
  read(1,*) L, alpha
  read(1,*)
  read(1,*) Nstep, delta_tau
  read(1,*)
  read(1,*) Nsample
  close(1)
 
  write(6,'(a)')        '================ calculation condition ================='
  write(6,'(a,I8)')     '  L (system size)                    =', L
  write(6,'(a,I8)')     '  alpha (hidden variable density)    =', alpha
  write(6,'(a,I8)')     '  total number of optimization steps =', Nstep
  write(6,'(a,F14.10)') '  delta tau                          =', delta_tau 
  write(6,'(a,I8)')     '  number of Monte Carlo samples      =', Nsample
  write(6,'(a)')        '========================================================'
  write(6,*)

  N  = L          ! number of visible units 
  M  = alpha * N  ! number of hidden units
  Nv = alpha * N  ! number of variational parameters
  Nwarmup = Nsample / 20
  ! 
  ! allocate arrays 
  !  
  allocate( Wirr(0:N-1,alpha) ); Wirr = 0d0
  allocate( x (L)             ); x  = 0  
  allocate( xp(L)             ); xp = 0  
  allocate( theta (M)         ); theta  = 0d0 
  allocate( thetap(M)         ); thetap = 0d0 
  allocate( gvec(Nv)          ); gvec = 0d0
  allocate( Ovec(Nv)          ); Ovec = 0d0
  allocate( Ovec_loc_x(Nv)    ); Ovec_loc_x = 0d0
  allocate( Smat(Nv,Nv)       ); Smat = 0d0
  ! 
  ! initialize W (variational parameters)   
  ! put random numbers between -0.01 and 0.01
  !
  open(unit=1,file='initial_W.txt',status='unknown')
  write(1,'(a)') '# index for independent neurons (f), index for interaction (iw), irreducible interaction Wirr(iw,f)'
  do f  = 1, alpha     
    do iw = 0, N-1 
      Wirr(iw,f) = 0.02d0*(grnd()-0.5d0)
      write(1,'(2I5,F22.15)') f, iw, Wirr(iw,f)
    end do ! iw 
    write(1,*)
  end do ! f 
  close(1)  
  ! 
  ! start optimization  
  ! 
  open(unit=10,file='Energy_vs_Iteration.txt',status='unknown')
  do iteration = 1, Nstep 
    ! 
    ! initialize configuration 
    ! (only for the 1st iteration, 
    !  from the 2nd iteration, we reuse configurations from the previous iteration)
    ! 
    if( iteration == 1 ) then
      do i = 1, N/2
        x(2*i-1) =  1
        x(2*i)   = -1
      end do ! i
    end if
    call calc_theta(x,theta)
    call calc_amplitude_RBM(M,theta,psi_x)
    ! 
    ! measure energy and calculate derivative of energy with respect to W
    !
    ! <H> = \sum_x p_x Eloc_x
    !
    !   p_x = | psi(x) |**2
    !   Eloc_x = \sum_x'  <x|H|x'> * ( psi(x')/psi(x) )
    !
    !  
    ! gvec(k) = 2 <H*O_k> - 2 <H> <O_k>   
    !   (for details, see e.g., Tahara and Imada, J. Phys. Soc. Jpn. 77, 114701 (2008)) 
    ! 
    !   <H*O_k> = \sum_x p_x Eloc_x Ovec_loc_x(k) 
    !   Ovec(k) = <O_k> = \sum_x p_x Ovec_loc_x(k) 
    !   
    !
    ! Smat(k1,k2) = <O_k1*O_k2> - <O_k1><O_k2>
    !
    gvec = 0d0
    Ovec = 0d0
    E    = 0d0
    Smat = 0d0
    do isample = 1, Nsample+Nwarmup
      ! 
      ! sample updates
      ! 
      do iupdate = 1, N
        call spin_flip_candidate(x,i1,i2)
        xp(:) = x(:)
        xp(i1) = -xp(i1)
        xp(i2) = -xp(i2)
        call calc_theta(xp,thetap)
        call calc_amplitude_RBM(M,thetap,psi_xp)
        if( grnd() > (psi_xp/psi_x)**2 ) cycle ! reject
        ! 
        ! accept => update configuration
        ! 
        x(:) = xp(:)
        theta(:) = thetap(:)
        psi_x = psi_xp
      end do ! iupdate
      if( isample <= Nwarmup ) cycle
      ! 
      ! calculate Eloc(x)
      ! H = J ( SzSz - SxSx - SySy ) after the gauge transformation
      !
      Eloc_x = 0d0
      do i1 = 1, N
        i2 = i1 + 1
        if( i1 == N ) i2 = 1
        !
        ! SzSz contribution 
        ! 
        Eloc_x = Eloc_x + x(i1)*x(i2) 
        !
        ! ( SxSx + SySy ) = 2 * ( S+S- + S-S+ ) contribution 
        ! 
        if( x(i1) /= x(i2) ) then 
          xp(:) = x(:)
          xp(i1) = -xp(i1)
          xp(i2) = -xp(i2)
          call calc_theta(xp,thetap)
          call calc_amplitude_RBM(M,thetap,psi_xp)
          Eloc_x = Eloc_x - 2d0*psi_xp/psi_x 
        end if
         
      end do ! i1  
      !
      ! calculate Ovec_loc_x 
      !
      k = 0 
      do f  = 1, alpha
      do iw = 0, N-1  

        k = k + 1
        Ovec_loc_x(k) = 0d0

        do jj = 1, N 
          j = (f-1)*N+jj
          i = jj + iw
          if( i > N ) i = i - N 
          Ovec_loc_x(k) = Ovec_loc_x(k) + tanh(theta(j))*dble(x(i))
        end do ! jj 

      end do ! iw
      end do ! f 
      if( k /= Nv ) stop 'k /= Nv'
      !
      ! summation for E, gvec, and Ovec 
      !
      E = E + Eloc_x  
      do k = 1, Nv
        gvec(k) = gvec(k) + 2d0 * Eloc_x * Ovec_loc_x(k)  
        Ovec(k) = Ovec(k) + Ovec_loc_x(k)  
      end do ! k 
      !
      ! summation for Smat
      !
      do k2 = 1, Nv
      do k1 = 1, Nv
        Smat(k1,k2) = Smat(k1,k2) + Ovec_loc_x(k1) * Ovec_loc_x(k2)  
      end do ! k1 
      end do ! k2 

    end do ! isample
    ! 
    ! By taking Monte Carlo average, we get  
    !  gvec(k) = 2 <H*O_k> 
    !  Ovec(k) = <O_k>  
    !  E = <H>  
    !  Smat(k1,k2) = <O_k1*O_k2>
    !  
    gvec(:) = gvec(:) / dble(Nsample)
    Ovec(:) = Ovec(:) / dble(Nsample)
    E = E / dble(Nsample) 
    Smat(:,:) = Smat(:,:) / dble(Nsample)
    ! 
    ! gvec(k) = 2 <H*O_k> - 2 <H> <O_k>  
    !  
    do k = 1, Nv
      gvec(k) = gvec(k) - 2d0 * E * Ovec(k)
    end do ! k
    ! 
    ! Smat(k1,k2) = <O_k1*O_k2> - <O_k1><O_k2>
    !  
    Sdiag_av = 0d0
    do k2 = 1, Nv
    do k1 = 1, Nv
      Smat(k1,k2) = Smat(k1,k2) - Ovec(k1) * Ovec(k2)  
      if( k1 == k2 ) then 
        if( Smat(k1,k2) <= 0d0 ) stop 'diagonal element of S matrix is not positive'
        Sdiag_av = Sdiag_av + Smat(k1,k2)
      end if 
    end do ! k1 
    end do ! k2 
    Sdiag_av = Sdiag_av / dble(Nv)
    !
    ! stabilization factor 
    !   For details, see Appendix C in Nomura et al., Phys. Rev. B 96, 205152 (2017) 
    !   Here, we add small uniform constant to diagonal elements of S matrix.
    !   A more aggressive optimization sometimes helps to lower the energy.
    !
    do k1 = 1, Nv
      Smat(k1,k1) = Smat(k1,k1) + Sdiag_av*0.0001d0 
    end do ! k1
    ! 
    ! update W (variational parameters) 
    ! 
    call solve_linear_equation(Nv,1,Smat,gvec) 
    delW_max = 0d0 
    do k = 1, Nv
      if( abs(delta_tau*gvec(k)) > delW_max ) delW_max = abs(delta_tau*gvec(k))
    end do ! k 
    k = 0 
    do f  = 1, alpha     
    do iw = 0, N-1 
      k = k + 1
      if( delW_max > 0.03d0 ) then  
        ! avoid large change in variational parameters to stabilize optimization
        Wirr(iw,f) = Wirr(iw,f) - delta_tau*gvec(k) * (0.03d0/delW_max) 
      else 
        Wirr(iw,f) = Wirr(iw,f) - delta_tau*gvec(k) 
      end if
    end do ! iw
    end do ! f
    ! 
    ! write information 
    ! 
    write(6,'(a,I6)')     'Iteration :', iteration
    write(6,'(a,F20.10)') '   total energy        :', E
    write(6,'(a,F20.10)') '   maximum change in W :', min(delW_max,0.03d0) 
    write(6,'(a)')        '   W parameters:'
    do f = 1, alpha 
      write(6,'(5x,20F9.5)')  (Wirr(iw,f), iw = 0, N-1)
    end do ! f

    write(10,*) iteration, E

  end do ! iteration
  close(10)
  ! 
  ! write optimized variational parameters
  ! 
  open(unit=1,file='optimized_W.txt',status='unknown')
  write(1,'(a)') '# index for independent neurons (f), index for interaction (iw), irreducible interaction Wirr(iw,f)'
  do f  = 1, alpha     
    do iw = 0, N-1 
      write(1,'(2I5,F22.15)') f, iw, Wirr(iw,f)
    end do ! iw
    write(1,*) 
  end do ! f 
  close(1)  

  deallocate(x,xp,Wirr,theta,thetap,gvec,Smat,Ovec,Ovec_loc_x)
end program 
