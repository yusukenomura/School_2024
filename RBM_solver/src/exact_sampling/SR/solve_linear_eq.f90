!
! solving A*x = b  
!
subroutine solve_linear_equation(N,N_rhs,A,B) 
  implicit none 
  integer , intent(in)    :: N 
  integer , intent(in)    :: N_rhs
  real(8) , intent(inout) :: A(N,N) 
  real(8) , intent(inout) :: B(N,N_rhs) 
  integer :: info
  integer :: i, j 
  real(8) :: r1

  do j = 1, N 
  do i = 1, N 
    r1 = A(i,j) - A(j,i)
    if( abs(r1) > 1.0d-5 ) write(6,*) & 
        'Warning: matrix A is not symmetric (solve_linear_equation)'
  end do ! i 
  end do ! j  

  info = 0
  call dposv('U',N,N_rhs,A,N,B,N,info)  

  if(info /= 0) then
    write(6,*) 'info (subrouitine solve_linear_equation):' , info
    stop
  end if 

  return 
end subroutine 
