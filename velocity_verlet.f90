module velocity_verlet
  ! module to solve Newton mechanics for a N-body planetary system
  ! using velocity verlet integration
  implicit none
  
  public :: verlet_n
  
  private :: compute_new_position, compute_half_step_velocity, &
       compute_new_velocity
  
contains
  
  subroutine verlet_n(bodies,nbodies,n,dt)
    ! solver for a N-body system
    use body
    implicit none
    integer, intent(in) :: nbodies
    type(planetary_body), dimension(nbodies), intent(inout) :: bodies
    integer, intent(in) :: n
    real(8), intent(in) :: dt
    real(8), dimension(nbodies) :: m
    real(8), dimension(nbodies,3) :: vh, r, newforces
    real(8), dimension(3) :: pos, vel, f
    integer :: z
    ! half step velocities and new positions
    do z=1,nbodies
       call get_body_fields(bodies(z),n,pos,vel,f,m(z))
       vh(z,:) = compute_half_step_velocity(vel,f,m(z),dt)
       r(z,:) = compute_new_position(pos,vh(z,:),dt)
    end do
    ! new forces and velocities
    newforces = calc_forces(m,r,nbodies)
    do z=1,nbodies
       bodies(z)%f(n+1,:) = newforces(z,:)
       bodies(z)%vel(n+1,:) = compute_new_velocity(bodies(z)%f(n+1,:),m(z),vh(z,:),dt)
       bodies(z)%pos(n+1,:) = r(z,:)
    end do        
  end subroutine verlet_n
  
  function compute_half_step_velocity(v0,f0,m,dt) result(v_half)
    ! Velocity at half time step from initial
    ! velocity and force
    implicit none
    real(8), dimension(3), intent(in) :: v0, f0
    real(8), intent(in) :: m, dt
    real(8), dimension(3) :: v_half
    v_half = v0 + 0.5*f0/m*dt
  end function compute_half_step_velocity

  function compute_new_position(r0,v_half,dt) result(r1)
    ! New position from initial position
    ! and half step velocity
    implicit none
    real(8), dimension(3), intent(in) :: r0, v_half
    real(8), intent(in) :: dt
    real(8), dimension(3) :: r1
    r1 = r0 + v_half*dt
  end function compute_new_position
  
  function compute_new_velocity(f1,m,v_half,dt) result(v1)
    ! Velocity at full step
    real(8), dimension(3), intent(in) :: f1, v_half
    real(8), intent(in) :: m, dt    
    real(8), dimension(3) :: v1
    v1 = v_half + 0.5*f1/m*dt
  end function compute_new_velocity  

end module velocity_verlet
