!
! compute L-BFGS direction
!
!  ndim        : dimension of the problem 
!  mhist       : number of historical iterations
!  iter        : outer iteration number, should start from 1.
!  x0          : current x
!  g0          : gradient at x
!  x(:,1:m+1)  : working array, store m+1 previous x
!  g(:,1:m+1)  : working array, store m+1 previous g

!
! output:
!  q: L-BFGS direction
!  x: updated on exit
!  g: updated on exit
!
! based on
!   https://github.com/GuipengLi/optLBFGS/blob/master/optLBFGS.m
!
!   See Alg. 3 in "A Stochastic Quasi-Newton Method for Online Convex Optimization", N.N. Schraudolph et al.
!     Proceedings of Machine Learning Research, 2:436-443, (2007)
!     http://proceedings.mlr.press/v2/schraudolph07a.html
!
! created on 1/21/2019 by Chen Huang
!
subroutine lbfgs_dir(ndim,mhist,iter,g0,x0,g,x,q)
  implicit none
  integer, intent(in) :: ndim, mhist, iter
  real(8)             :: g(ndim,mhist+1), &
                         x(ndim,mhist+1), &
                         g0(ndim), x0(ndim)

  real(8), intent(out) :: q(ndim)

  ! local vars
  integer :: i, j, m
  real(8) :: s(ndim), y(ndim), rho, alpha(mhist), beta


  ! >>>> FUNCTION BEGINS <<<<

  ! following Wikipedia
  ! https://en.wikipedia.org/wiki/Limited-memory_BFGS#Algorithm

  ! update historical x and g
  ! 1,2,...,m are historical x and g
  do i=1,mhist
    x(:,i) = x(:,i+1)
    g(:,i) = g(:,i+1)
  enddo
  ! m+1 is the current x and g
  x(:,mhist+1) = x0
  g(:,mhist+1) = g0


  if (mhist >= iter-1) then
    m = iter-1
  else
    m = mhist
  endif


  if (m==0) then
    q = -g0
    return
  endif

  q = -g0

  ! first loop
  do i = 1, m
    s = x(:,mhist-i+2)-x(:,mhist-i+1)  ! s_i
    y = g(:,mhist-i+2)-g(:,mhist-i+1)  ! y_i
    rho = 1.d0/sum(y*s)      ! define rho
    alpha(i) = rho*sum(s*q)  ! define alpha
    q = q - alpha(i)*y
  enddo

  s = x0 - x(:,mhist)
  y = g0 - g(:,mhist)
  q = sum(s*y)/sum(y*y)*q

  ! second loop
  do i = m,1,-1
    s = x(:,mhist-i+2)-x(:,mhist-i+1)  ! s_i
    y = g(:,mhist-i+2)-g(:,mhist-i+1)  ! y_i
    rho = 1.d0/sum(y*s)      ! define rho
    beta = rho*sum(y*q)      ! define beta
    q = q + s*(alpha(i)-beta)
  enddo

end subroutine lbfgs_dir
