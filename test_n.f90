program test_n
  ! program to simulate planetary system with N objects
  use body
  use velocity_verlet
  implicit none

  integer :: nsteps
  real(8) :: duration, step
  real(8), dimension(:,:), allocatable :: vel, pos
  real(8), dimension(:), allocatable :: mass 
  type(planetary_body), dimension(:), allocatable :: bodies
  integer :: n, nbodies
  character(len=1024) :: input_file
  
  input_file = 'input.txt'

  ! read input values
  call read_input_parameters(input_file,nbodies,duration,step,mass,vel,pos)

  ! allocate planetary objects
  allocate(bodies(nbodies))

  ! number of steps in this simulation
  nsteps = ceiling(duration*60.0*60.0*24.0/step)

  ! init objects
  do n=1,nbodies
     call init_body(bodies(n),pos(n,:),vel(n,:),mass(n),nsteps) 
  end do

  ! initial forces
  call init_n_body_system(bodies,nbodies)

  ! run simulation 
  do n=1,nsteps
     call verlet_n(bodies,nbodies,n,step)
  end do
  
  ! write output
  call write_output(bodies)

contains 

  subroutine read_input_parameters(fname,nbodies,duration,step,mass,vel,pos)    
    ! Read input values from a file
    implicit none
    character(len=1024), intent(in) :: fname
    integer, intent(out) :: nbodies
    real(8), intent(out) :: duration, step
    real(8), dimension(:,:), allocatable, intent(out) :: vel, pos
    real(8), dimension(:), allocatable, intent(out) :: mass 
    open(unit=10,file=fname,status='old')
    read(10,'(/,I6,///,F10.3,///,F10.3,//)') nbodies, duration, step
    allocate(mass(nbodies), vel(nbodies,3), pos(nbodies,3))
    do n=1,nbodies
       read(10,*) mass(n)
    end do
    read(10,'(/)')
    do n=1,nbodies
       read(10,*) pos(n,:)
    end do
    read(10,'(/)')
    do n=1,nbodies
       read(10,*) vel(n,:)
    end do    
    close(10)
  end subroutine read_input_parameters

  subroutine write_output(bodies)
    ! writes positions for each time steps 
    ! for all objects
    implicit none
    type(planetary_body), dimension(:), intent(in) :: bodies
    integer :: n, m
    character(len=1024) :: ofile
    ! write positions
    do n=1, size(bodies,1)
       write(ofile,'(I1)') n
       ofile = 'results/pos_' // trim(ofile) // '.dat'
       open(unit=10,file=ofile)
       do m=1,size(bodies(n)%pos,1)
          write(10,*) bodies(n)%pos(m,:)
       end do
       close(10)
    end do
    ! write number of steps in a separate file
    open(unit=10,file='results/n_steps.dat')
    write(10,*) size(bodies(1)%pos,1)
    close(10)    
  end subroutine write_output

end program test_n
