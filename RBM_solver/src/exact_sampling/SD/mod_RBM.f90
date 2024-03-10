module mod_RBM
  integer :: N                       ! number of visible units
  integer :: M                       ! number of hidden units
  integer :: alpha                   ! hidden variable density = M/N 
  real(8) , allocatable :: Wirr(:,:) ! irreducible part of RBM interaction parameters
end module
