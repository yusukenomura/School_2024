!
!  Solving 1D antiferromagnetic Heisenberg model by RBM 
!    using steepest descent (SD) for optimization  
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

  integer :: D         ! dimension of Hamiltonian matrix
  integer :: L         ! system size
  integer :: Nv        ! number of variational parameters 
  integer :: Nstep     ! total number of optimization steps
  real(8) :: delta_tau ! controls the magnitude of parameter update 

  integer :: f         ! index for independent neurons 
  integer :: i         ! index for visible units
  integer :: j         ! index for hidden units
  integer :: k         ! index for variational parameters 
  integer :: jj        ! subindex for hidden units
  integer :: ix        ! index for spin configurations
  integer :: iw        ! index for W interaction 
  integer :: iteration ! index for optimization iterations

  real(8) :: psi_x     ! psi(x)  = <x |psi>
  real(8) :: psi_xp    ! psi(x') = <x'|psi>
  real(8) :: norm      ! <psi|psi>
  real(8) :: p_x       ! p(x) = | psi(x) |**2
  real(8) :: Eloc_x    ! Eloc(x) = \sum_x'  <x|H|x'> * ( psi(x')/psi(x) )
  real(8) :: E         ! total energy 
  real(8) :: delW_max  ! maximum change in W at each optimization step

  integer , allocatable :: xlist(:,:)    ! list of spin configurations
  integer , allocatable :: x (:)         ! |x>  = |sigma_1 , sigma_2 , ..., sigma_N >
  integer , allocatable :: xp(:)         ! |x'> = |sigma'_1, sigma'_2, ..., sigma'_N>
  real(8) , allocatable :: theta (:)     ! theta (j) = \sum_i  W_ij * \sigma_i
  real(8) , allocatable :: thetap(:)     ! thetap(j) = \sum_i  W_ij * \sigma'_i 
  real(8) , allocatable :: gvec(:)       ! derivative of energy with respect to variational parameters 


  real(8) :: r1 
  integer :: i1, i2
  integer :: seedsize
  integer , allocatable :: seed(:)
  real(8) , allocatable :: Ovec(:)        
  real(8) , allocatable :: Ovec_loc_x(:)  
  real(8) , allocatable :: psi(:)
  ! 
  ! set seed for random number generator  
  !
  call random_seed(size=seedsize)
  allocate( seed(seedsize) )
  do i = 1, seedsize 
    call system_clock(count=seed(i))
  end do ! i
  call random_seed(put=seed(:))
  ! 
  ! read inputs
  ! 
  open(unit=1,file='RBM.input',status='old')
  read(1,*)
  read(1,*) alpha, Nstep, delta_tau
  close(1)

  open(unit=1,file='spin_configurations.txt',status='old')
  read(1,*)
  read(1,*) L, D
  read(1,*)
  allocate( xlist(L,D) ); xlist = 0  
  do ix = 1, D
    read(1,*) xlist(:,ix)    
  end do ! ix 
  close(1) 
 
  write(6,'(a)')        '================ calculation condition ================='
  write(6,'(a,I8)')     '  L (system size)                    =', L
  write(6,'(a,I8)')     '  alpha (hidden variable density)    =', alpha
  write(6,'(a,I8)')     '  total number of optimization steps =', Nstep
  write(6,'(a,F14.10)') '  delta tau                          =', delta_tau 
  write(6,'(a)')        '========================================================'
  write(6,*)

  N  = L          ! number of visible units 
  M  = alpha * N  ! number of hidden units
  Nv = alpha * N  ! number of variational parameters
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
  ! 
  ! initialize W (variational parameters)   
  ! put random numbers between -0.01 and 0.01
  !
  open(unit=1,file='initial_W.txt',status='unknown')
  write(1,'(a)') '# index for independent neurons (f), index for interaction (iw), irreducible interaction Wirr(iw,f)'
  do f  = 1, alpha     
    do iw = 0, N-1 
      call random_number(r1)
      Wirr(iw,f) = 0.02d0*(r1-0.5d0)
      write(1,'(2I5,F22.15)') f, iw, Wirr(iw,f)
    end do ! iw 
    write(1,*)
  end do ! f 
  close(1)  
  ! 
  ! initial wave function
  ! 
  allocate( psi(D) ); psi = 0d0
  norm = 0d0
  do ix = 1, D
    x(:) = xlist(:,ix)
    call calc_theta(x,theta)
    call calc_amplitude_RBM(M,theta,psi(ix))
    norm = norm + psi(ix)*psi(ix)
  end do ! ix 

  open(unit=1,file='initial_wf.txt',status='unknown')
  write(1,'(a)') '# initial wave function psi(x), order of x follows that in spin_configurations.txt'
  do ix = 1, D
    write(1,'(F25.15)') psi(ix) / dsqrt(norm) 
  end do ! ix 
  deallocate(psi)
  close(1)
  ! 
  ! start optimization  
  ! 
  open(unit=10,file='Energy_vs_Iteration.txt',status='unknown')
  do iteration = 1, Nstep 
    ! 
    ! measure energy and calculate derivative of energy with respect to W
    ! here, we do not use Monte Carlo, instead we take brute-force sum over x 
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
    norm = 0d0 
    gvec = 0d0
    Ovec = 0d0
    E    = 0d0
    do ix = 1, D
      ! 
      ! calculate wave function psi(x)
      ! 
      x(:) = xlist(:,ix)
      call calc_theta(x,theta)
      call calc_amplitude_RBM(M,theta,psi_x)
      p_x  = psi_x*psi_x
      norm = norm + psi_x*psi_x 
      ! 
      ! calculate Eloc(x) and take sum over x for E
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
      E = E + p_x * Eloc_x  
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
      ! summation over x for gvec and Ovec 
      !
      do k = 1, Nv
        gvec(k) = gvec(k) + 2d0 * p_x * Eloc_x * Ovec_loc_x(k)  
        Ovec(k) = Ovec(k) + p_x * Ovec_loc_x(k)  
      end do ! k 

    end do ! ix
    ! 
    ! Dividing by <psi|psi> (=norm), we get  
    !  gvec(k) = 2 <H*O_k> 
    !  Ovec(k) = <O_k>  
    !  E = <H>  
    !  
    gvec(:) = gvec(:) / norm
    Ovec(:) = Ovec(:) / norm
    E = E / norm 
    ! 
    ! gvec(k) = 2 <H*O_k> - 2 <H> <O_k>  
    !  
    do k = 1, Nv
      gvec(k) = gvec(k) - 2d0 * E * Ovec(k)
    end do ! k
    ! 
    ! update W (variational parameters) 
    ! 
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
  ! 
  ! optimized wave function
  ! 
  allocate( psi(D) ); psi = 0d0
  norm = 0d0
  do ix = 1, D
    x(:) = xlist(:,ix)
    call calc_theta(x,theta)
    call calc_amplitude_RBM(M,theta,psi(ix))
    norm = norm + psi(ix)*psi(ix)
  end do ! ix 

  open(unit=1,file='optimized_wf.txt',status='unknown')
  write(1,'(a)') '# optimized wave function psi(x), order of x follows that in spin_configurations.txt'
  do ix = 1, D
    write(1,'(F25.15)') psi(ix) / dsqrt(norm) 
  end do ! ix 
  deallocate(psi)
  close(1)

  deallocate(seed,xlist,x,xp,Wirr,theta,thetap,gvec,Ovec,Ovec_loc_x)
end program 
