module body
  ! - Module to create and initialize "planetary_body" type objects.
  ! - These objects are described by position, velocity, force, and mass.
  ! - Position, velocity and force are 3D-vectors. Their values are
  !   saved for each time step. 
  implicit none
  
  type planetary_body
     real(8), dimension(:,:), allocatable :: pos ! position (x,y,z) [m]
     real(8), dimension(:,:), allocatable :: vel ! velocity (dx,dy,dz) [m/s]
     real(8), dimension(:,:), allocatable :: f   ! force (dx,dy,dz) [N]
     real(8) :: mass                             ! mass [kg]
  end type planetary_body
  
contains

  subroutine create_body(body,n)
    ! create body type structure
    ! with n-steps of information + (inital state)
    implicit none
    integer, intent(in) :: n
    type(planetary_body), intent(inout) :: body
    if (allocated(body%pos)) deallocate(body%pos)
    if (allocated(body%vel)) deallocate(body%vel)
    if (allocated(body%f)) deallocate(body%f)
    allocate(body%pos(n+1,3))
    allocate(body%vel(n+1,3))
    allocate(body%f(n+1,3))
  end subroutine create_body

  subroutine init_body(body,pos,vel,mass,nsteps)
    ! initial state of the body (with mass, position, velocity)
    implicit none
    type(planetary_body), intent(inout) :: body
    real(8), dimension(3), intent(in) :: pos, vel
    real(8), intent(in) :: mass
    integer, intent(in) :: nsteps
    call create_body(body,nsteps)
    body%mass = mass
    body%pos(1,:) = pos
    body%vel(1,:) = vel
  end subroutine init_body

  subroutine get_body_fields(body,n,pos,vel,f,mass)
    ! returns position, velocity, force and mass of
    ! a object at the time step n
    implicit none
    integer, intent(in) :: n
    type(planetary_body), intent(inout) :: body
    real(8), dimension(3), intent(out) :: pos, vel, f
    real(8), intent(out) :: mass
    pos = body%pos(n,:)
    vel = body%vel(n,:)
    f = body%f(n,:)
    mass = body%mass
  end subroutine get_body_fields

  subroutine init_n_body_system(bodies,nbodies)
    ! computes initial forces of a N body system
    ! using Newton's universal law of gravity
    use math,only : newton
    implicit none
    integer, intent(in) :: nbodies
    type(planetary_body), dimension(nbodies), intent(inout) :: bodies
    real(8), dimension(nbodies,3) :: pos
    real(8), dimension(nbodies) :: mass
    integer :: n
    real(8), dimension(nbodies,3) :: f
    ! initial positions and masses
    do n=1, nbodies
       pos(n,:) = bodies(n)%pos(1,:)
       mass(n) = bodies(n)%mass
    end do
    ! initial forces
    f = calc_forces(mass,pos,nbodies)
    do n=1, nbodies
       bodies(n)%f(1,:) = f(n,:)
    end do
  end subroutine init_n_body_system

  function calc_forces(mass,pos,nbodies) result(f)
    ! calculate forces between N-bodies
    use math, only : newton
    implicit none
    integer, intent(in) :: nbodies
    real(8), dimension(nbodies), intent(in) :: mass
    real(8), dimension(nbodies,3), intent(in) :: pos
    real(8), dimension(nbodies,3) :: f
    integer :: n, m, i1, i2
    f = 0.0
    do n=1, nbodies
       do m=1, nbodies-1
          i1 = n
          i2 = m
          if (m>=n) i2 = i2+1
          f(n,:) = f(n,:) + newton(mass(i2),mass(i1),pos(i2,:),pos(i1,:))
       end do
    end do
  end function calc_forces
  
end module body
