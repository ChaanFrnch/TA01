!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!< Defines all variables needed to parametrize a simulation run.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module HydroParameters

  use HydroPrecision
  use HydroConstants
  use Partitioner
  use mpi

  !! run parameters
  integer(int_kind ) :: nStepmax !< maximun number of time steps.
  real   (fp_kind)   :: tEnd     !< end of simulation time.
  integer(int_kind ) :: nOutput  !< number of time steps between 2 consecutive outputs.

  !! geometry parameters
  integer(int_kind) :: nx=0     !< logical size along X (without ghost cells).
  integer(int_kind) :: ny=0     !< logical size along Y (without ghost cells).
  integer(int_kind), parameter :: ghostWidth=2 
  integer(int_kind) :: imin=0   !< index minimum at X border
  integer(int_kind) :: imax=0   !< index maximum at X border
  integer(int_kind) :: jmin=0   !< index minimum at Y border
  integer(int_kind) :: jmax=0   !< index maximum at Y border

  integer(int_kind) :: isize_tot=0  !< total size (in cell unit) along X direction with ghosts.
  integer(int_kind) :: jsize_tot=0  !< total size (in cell unit) along Y direction with ghosts.
  integer(int_kind) :: isize=0  !< local size (in cell unit) along X direction with ghosts.
  integer(int_kind) :: jsize=0  !< local size (in cell unit) along Y direction with ghosts.

  real(fp_kind)   :: xmin=0.0
  real(fp_kind)   :: xmax=1.0
  real(fp_kind)   :: ymin=0.0
  real(fp_kind)   :: ymax=1.0
  real(fp_kind)   :: dx       !< x resolution
  real(fp_kind)   :: dy       !< y resolution

  integer(int_kind) :: boundary_type_xmin=1
  integer(int_kind) :: boundary_type_xmax=1
  integer(int_kind) :: boundary_type_ymin=1
  integer(int_kind) :: boundary_type_ymax=1

  !! IO parameters
  logical :: ioVTK=.true.   !< enable VTK  output file format (using VTI).
  logical :: ioHDF5=.false. !< enable HDF5 output file format.

  !! hydro (numerical scheme) parameters
  real   (fp_kind) :: gamma0        !< specific heat capacity ratio (adiabatic index)
  real   (fp_kind) :: gamma6
  real   (fp_kind) :: cfl           !< Courant-Friedrich-Lewy parameter.
  real   (fp_kind) :: slope_type    !< type of slope computation (2 for second order scheme).
  integer(int_kind):: iorder        !< 
  real   (fp_kind) :: smallr        !< small density cut-off
  real   (fp_kind) :: smallc        !< small speed of sound cut-off
  real   (fp_kind) :: smallp        !< small pressure cut-off
  real   (fp_kind) :: smallpp       !< smallp times smallr
  integer(int_kind) :: niter_riemann=10 !< number of iteration usd in quasi-exact riemann solver
  integer(int_kind) :: riemannSolverType=0
  character(LEN=20) :: problem='implode'


  !! other parameters
  integer(int_kind)  ::  nbVar=4  !< number of fields in simulation (density, energy, vx, vy)
  integer(int_kind)  :: implementationVersion=0 !< triggers which implementation to use (currently 2 versions)
  integer :: coord_x, coord_y 
  integer :: size_x_max, size_y_max

  contains

    subroutine initHydroParameters(nbTask, myRank, sizes_x, sizes_y)
      
      implicit none

      ! local variables
      integer, intent(in) :: nbTask, myRank
      !integer, intent(in) :: coord_x, coord_y
      integer(int_kind) :: narg
      character(LEN=80) :: inputFilename
      integer :: size_x
      !! partitioner parameters
      integer :: nbcores, mx, my
      integer,dimension(:), intent(out) :: sizes_x, sizes_y
   
      ! declare namelist
      namelist/run/tEnd,nStepmax,nOutput
      namelist/mesh/nx,ny,boundary_type_xmin,boundary_type_xmax,boundary_type_ymin,boundary_type_ymax
      namelist/hydro/gamma0,cfl,smallr,smallc,niter_riemann, &
           &         iorder,slope_type
      namelist/other/implementationVersion

      narg = COMMAND_ARGUMENT_COUNT()
      if (narg .ne. 1) then
         write(*,*)'You should provide an input parameter file (with a namelist)'
         write(*,*)'Usage:'
         write(*,*)'   ./euler2d param.nml'
         stop
      end if

      ! retrieve input parameter file name
      call GET_COMMAND_ARGUMENT(1,inputFilename)

      ! read namelist
      open(1,file=inputFilename)
      read(1,NML=run)
      read(1,NML=mesh)
      read(1,NML=hydro)
      read(1,NML=other)
      close(1)

      ! set other parameters
      imin = 1
      jmin = 1
      imax = nx+2*ghostWidth
      jmax = ny+2*ghostWidth

      isize_tot = imax - imin + 1
      jsize_tot = jmax - jmin + 1

      nbcores = nbTask
      call partition (nbcores,nx,ny,mx,my,sizes_x,sizes_y)
      
      size_x_max = size(sizes_x)
      size_y_max = size(sizes_y)
      coord_x = modulo(myRank,size_x)
      coord_y = (myRank-coord_x)/size_x
      
      isize = sizes_x(coord_x+1) + 2*ghostwidth
      jsize = sizes_y(coord_y+1) + 2*ghostwidth
      write(*,*) 'I am proc ', myRank, 'and isize = ', isize, 'and jzise = ', jsize

      dx = (xmax - xmin) / nx
      dy = (ymax - ymin) / ny

      smallp  = smallc*smallc/gamma0
      smallpp = smallr*smallp
      gamma6  = (gamma0 + ONE_F)/(TWO_F * gamma0)

      ! check that given parameters are valid
      if ( (implementationVersion .ne. 0) .and. (implementationVersion .ne. 1)) then
         write(*,*) 'The implementation version parameter should 0 or 1 !!!'
         write(*,*) 'Check your namelist file, section OTHER'
         stop
      else
         write(*,*) 'Using implementation version', implementationVersion
      end if

    end subroutine initHydroParameters
    
    ! dump hydro parameters on stdout
    subroutine printHydroParameters()

      implicit none

      write(*,*) 'Simulation run parameters:'
      write(*,*) 'nx : ', nx
      write(*,*) 'ny : ', ny
      write(*,*) 'dx : ', dx
      write(*,*) 'dy : ', dy
      write(*,*) 'imin : ', imin, 'imax : ', imax
      write(*,*) 'jmin : ', jmin, 'jmax : ', jmax      
      write(*,*) 'nStepmax,tEnd,nOutput : ',nStepmax,tEnd,nOutput
      write(*,*) 'gamma0,cfl : ',gamma0,cfl
      write(*,*) 'smallr,smallc,niter_riemann :',smallr,smallc,niter_riemann
      write(*,*) 'iorder,slope_type :',iorder,slope_type

    end subroutine printHydroParameters

end module HydroParameters
