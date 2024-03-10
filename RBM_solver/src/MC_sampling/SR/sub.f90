!
! Subroutine for calculating theta 
!   see Eq.(S9) in Supplementary Materials, Carleo and Troyer, Science (2017)
!   theta(j) = \sum_i  W_ij * \sigma_i 
! 
subroutine calc_theta(sigma,theta)
  use mod_RBM , only : alpha, N, M, Wirr
  implicit none
  integer , intent(in)  :: sigma(N)
  real(8) , intent(out) :: theta(M)
  real(8) :: W(N,M)
  integer :: f, j, i, jj, iw

  j = 0 
  do f  = 1, alpha  ! loop for independent neuron
  do jj = 1, N      ! we have N copies of hidden units because of translational symmetry 

    j = j + 1   
    theta(j) = 0d0

    do i = 1, N
      ! 
      ! set W_ij using translational symmetry  
      !
      iw = i - jj  
      if( iw < 0 ) iw = iw + N  
      W(i,j) = Wirr(iw,f) 
      ! 
      ! calculate theta 
      !  
      theta(j) = theta(j) + W(i,j) * dble(sigma(i))
    end do ! i 

  end do ! jj 
  end do ! f 

  return 
end subroutine 
!
! Subroutine for calculating amplitude of RBM wave function
! 
subroutine calc_amplitude_RBM(M,theta,psi_x)
  implicit none
  integer , intent(in)  :: M        ! number of hidden units
  real(8) , intent(in)  :: theta(M) ! theta(j) = \sum_i  W_ij * \sigma_i
  real(8) , intent(out) :: psi_x    ! <x|psi>   (|x> = |sigma_1, sigma_2, ..., sigma_N>) 
  integer :: j
  ! 
  ! When bias terms are zero, wave function is given by 
  !   psi(x) = \prod_j 2d0 * cosh(theta(j)) 
  ! 
  psi_x = 1d0
  do j = 1, M
    psi_x = psi_x * cosh(theta(j)) ! neglect irrelevant factor of 2
  end do ! j

  return 
end subroutine 
! 
! x => xp (flip i1th and i2th spins)  
! 
subroutine spin_flip_candidate(x,i1,i2)
  use mod_RBM , only : N
  implicit none
  integer , intent(in)  :: x(N)
  integer , intent(out) :: i1, i2
  real(8) , external :: grnd
  logical :: lsuccess
  integer :: itmp

  lsuccess = .false. 
  do itmp = 1, 10000
    i1 = mod(int(dble(N)*grnd()),N) + 1
    i2 = mod(int(dble(N)*grnd()),N) + 1
    if( x(i1) /= x(i2) ) then
      lsuccess = .true.
      exit
    end if
  end do ! itmp

  if( .not. lsuccess ) stop 'subroutine spin_flip_candidate fails'

  return 
end subroutine 
