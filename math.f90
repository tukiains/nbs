module math
  ! Module with some mathematical functions
  implicit none

  public :: newton
  private :: norm, distance_between_points

contains
  
  function norm(r_in) result(norm_value)
    ! 2-norm of r
    implicit none
    real(8), dimension(3), intent(in) :: r_in
    real(8) :: norm_value
    norm_value = sqrt(sum((r_in)**2.0))
  end function norm
  
  function distance_between_points(r1,r2) result(r)
    ! distance between two points
    implicit none
    real(8), dimension(3), intent(in) :: r1, r2
    real(8) :: r
    r = norm(r2-r1)   
  end function distance_between_points
  
  function newton(m1,m2,r1,r2) result(f)
    ! Newton's universal law of gravity
    implicit none
    real(8), intent(in) :: m1, m2
    real(8), dimension(3), intent(in) :: r1, r2
    real(8), dimension(3) :: f, r12
    real(8) :: r
    real(8), parameter :: G = 6.67384e-11
    r = distance_between_points(r1,r2)
    r12 = (r2-r1)/r
    f = -G*m1*m2/r**2.0*r12
  end function newton
  
end module math
