!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!< Defines all variables needed to perform computations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module HydroRun

  use HydroPrecision
  use HydroConstants
  use HydroParameters
  use HydroUtils
  use Communication
  use mpi

  ! defines data arrays
  real(fp_kind), dimension(:,:,:), allocatable :: u,u2,u_tot !< conservative variables
  real(fp_kind), dimension(:,:,:), allocatable :: q    !< primitive variables (implementation version 1 only)
  real(fp_kind), dimension(:,:,:), allocatable :: qm_x, qm_y !< input to Riemann solvers (implementation version 1 only)
  real(fp_kind), dimension(:,:,:), allocatable :: qp_x, qp_y !< input to Riemann solvers (implementation version 1 only)
  integer(int_kind) :: nStep=0

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initHydroRun

    implicit none

    !integer, intent(in) :: myRank

    ! memory allocation)
    allocate( u (isize, jsize, nbVar) )
    allocate( u2(isize, jsize, nbVar) )
    allocate( u_tot (isize_tot, jsize_tot, nbVar) )

    if (implementationVersion .eq. 1) then
       allocate( q   (isize, jsize, nbVar) )
       allocate( qm_x(isize, jsize, nbVar) )
       allocate( qm_y(isize, jsize, nbVar) )
       allocate( qp_x(isize, jsize, nbVar) )
       allocate( qp_y(isize, jsize, nbVar) )
    end if

    ! initialize u (at t=0)
    select case (problem)
    case('implode') ! discontinuity line along the domain diagonal
       call init_implode(u)
    case default
       write(*,*) 'Unknown problem; default to implode'
       call init_implode(u)
    end select

  end subroutine initHydroRun

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cleanupHydroRun

    implicit none

    ! memory free
    deallocate(u,u2)
    deallocate(u_tot)

    if (implementationVersion .eq. 1) then
       deallocate(q,qm_x,qm_y,qp_x,qp_y) 
    end if

  end subroutine cleanupHydroRun

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute time step reduction
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_dt(dt, useU)
    implicit none

    ! dummy variables
    real   (fp_kind) , intent(out)  :: dt
    integer(int_kind), intent(in)   :: useU

    ! local variables
    real(fp_kind) :: invDt=0.0
    real(fp_kind) :: vx,vy,dt_min
    integer :: i,j, ierr
    real(fp_kind), dimension(nbVar) :: qLoc
    real(fp_kind)                   :: c

    ! for loop over inner region
    if (useU .eq. 0) then

       do j=ghostWidth+1,jsize-ghostWidth-1
          do i=ghostWidth+1,isize-ghostWidth-1

             call computePrimitives(u, i, j, c, qLoc)
             vx = c + abs(qLoc(IU))
             vy = c + abs(qLoc(IV))
             invDt = max(invDt, vx/dx + vy/dy)

          end do
       end do

    else

       do j=ghostWidth+1,jsize-ghostWidth-1
          do i=ghostWidth+1,isize-ghostWidth-1

             call computePrimitives(u2, i, j, c, qLoc)
             vx = c + abs(qLoc(IU))
             vy = c + abs(qLoc(IV))
             invDt = max(invDt, vx/dx + vy/dy)

          end do
       end do
    end if

    dt = cfl / invDt

    ! to determine the most constraining dt and return it to each processor
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(dt, dt_min, 1, MPI_REAL, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt_min, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    dt = dt_min
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  end subroutine compute_dt

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Wrapper to the actual computation routine
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine godunov_unsplit(dt)

    implicit none

    ! dummy variables
    real(fp_kind), intent(in) :: dt

    if ( modulo(nStep,2) .eq. 0 ) then
       call godunov_unsplit_cpu(u , u2, dt)
    else
       call godunov_unsplit_cpu(u2, u , dt)
    end if

  end subroutine godunov_unsplit

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Actual CPU computation of Godunov scheme
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine godunov_unsplit_cpu(data_in, data_out, dt)

    use Monitoring

    implicit none

    ! dummy variables
    real(fp_kind), dimension(isize, jsize, nbVar), intent(inout) :: data_in
    real(fp_kind), dimension(isize, jsize, nbVar), intent(inout) :: data_out
    real(fp_kind)                  , intent(in)    :: dt

    ! local variables
    integer :: i, j, ii, jj, iVar

    real(fp_kind) :: dtdx
    real(fp_kind) :: dtdy

    ! Local variables for trace computation
    ! we need to store qm/qp for current position i,j and i-1,j and i,j-1 
    ! that is 1+2=3 positions in total
    real(fp_kind), dimension(3, nbVar) :: qm_x_loc
    real(fp_kind), dimension(3, nbVar) :: qm_y_loc

    real(fp_kind), dimension(3, nbVar) :: qp_x_loc
    real(fp_kind), dimension(3, nbVar) :: qp_y_loc

    real(fp_kind), dimension(nbVar)    :: qLoc ! local primitive variables
    real(fp_kind), dimension(4, nbVar) :: qNeighbors

    real(fp_kind), dimension(2, nbVar) :: qm, qp
    real(fp_kind)                      :: c, cPlus, cMinus
    integer :: pos

    ! Local variables for Riemann problems solving
    real(fp_kind), dimension(nbVar) :: qleft
    real(fp_kind), dimension(nbVar) :: qright
    real(fp_kind), dimension(nbVar) :: qgdnv
    real(fp_kind), dimension(nbVar) :: flux_x
    real(fp_kind), dimension(nbVar) :: flux_y

    dtdx = dt / dx
    dtdy = dt / dy

    ! fill ghost cell in data_in
    call timerStart(boundaries_timer)
    call make_boundaries(data_in)
    call timerStop(boundaries_timer)

    ! copy data_in into data_out (not necessary)
    data_out = data_in

    ! start main computation
    call timerStart(godunov_timer)

    if (implementationVersion .eq. 0) then

       do j=ghostWidth+1, jsize-ghostWidth+1
          do i=ghostWidth+1, isize-ghostWidth+1

             ! compute qm, qp for the 1+2 positions
             do pos=1,3

                ii=i
                jj=j
                if (pos .eq. 2) ii = i-1
                if (pos .eq. 3) jj = j-1
                call computePrimitives(data_in, ii  , jj  , c     , qLoc)
                call computePrimitives(data_in, ii+1, jj  , cPlus , qNeighbors(1,:))
                call computePrimitives(data_in, ii-1, jj  , cMinus, qNeighbors(2,:))
                call computePrimitives(data_in, ii  , jj+1, cPlus , qNeighbors(3,:))
                call computePrimitives(data_in, ii  , jj-1, cMinus, qNeighbors(4,:))

                ! compute qm, qp
                call trace_unsplit_2d(qLoc, qNeighbors, dtdx, dtdy, qm, qp)

                ! store qm, qp
                do iVar=1,nbVar
                   qm_x_loc(pos,iVar) = qm(1,iVar)
                   qp_x_loc(pos,iVar) = qp(1,iVar)
                   qm_y_loc(pos,iVar) = qm(2,iVar)
                   qp_y_loc(pos,iVar) = qp(2,iVar)
                end do ! end do iVar

             end do !! end do pos

             ! Solve Riemann problem at X-interfaces and compute X-fluxes
             qleft(ID)   = qm_x_loc(2,ID)
             qleft(IP)   = qm_x_loc(2,IP)
             qleft(IU)   = qm_x_loc(2,IU)
             qleft(IV)   = qm_x_loc(2,IV)

             qright(ID)  = qp_x_loc(1,ID)
             qright(IP)  = qp_x_loc(1,IP)
             qright(IU)  = qp_x_loc(1,IU)
             qright(IV)  = qp_x_loc(1,IV)

             call riemann_2d(qleft,qright,qgdnv,flux_x)

             ! Solve Riemann problem at Y-interfaces and compute Y-fluxes
             qleft(ID)   = qm_y_loc(3,ID)
             qleft(IP)   = qm_y_loc(3,IP)
             qleft(IU)   = qm_y_loc(3,IV) ! watchout IU, IV permutation
             qleft(IV)   = qm_y_loc(3,IU) ! watchout IU, IV permutation

             qright(ID)  = qp_y_loc(1,ID)
             qright(IP)  = qp_y_loc(1,IP)
             qright(IU)  = qp_y_loc(1,IV) ! watchout IU, IV permutation
             qright(IV)  = qp_y_loc(1,IU) ! watchout IU, IV permutation

             call riemann_2d(qleft,qright,qgdnv,flux_y)

             !
             ! update hydro array
             !
             data_out(i-1,j  ,ID) = data_out(i-1,j  ,ID) - flux_x(ID)*dtdx
             data_out(i-1,j  ,IP) = data_out(i-1,j  ,IP) - flux_x(IP)*dtdx
             data_out(i-1,j  ,IU) = data_out(i-1,j  ,IU) - flux_x(IU)*dtdx
             data_out(i-1,j  ,IV) = data_out(i-1,j  ,IV) - flux_x(IV)*dtdx

             data_out(i  ,j  ,ID) = data_out(i  ,j  ,ID) + flux_x(ID)*dtdx
             data_out(i  ,j  ,IP) = data_out(i  ,j  ,IP) + flux_x(IP)*dtdx
             data_out(i  ,j  ,IU) = data_out(i  ,j  ,IU) + flux_x(IU)*dtdx
             data_out(i  ,j  ,IV) = data_out(i  ,j  ,IV) + flux_x(IV)*dtdx

             data_out(i  ,j-1,ID) = data_out(i  ,j-1,ID) - flux_y(ID)*dtdx
             data_out(i  ,j-1,IP) = data_out(i  ,j-1,IP) - flux_y(IP)*dtdx
             data_out(i  ,j-1,IU) = data_out(i  ,j-1,IU) - flux_y(IV)*dtdx ! watchout IU and IV swapped
             data_out(i  ,j-1,IV) = data_out(i  ,j-1,IV) - flux_y(IU)*dtdx ! watchout IU and IV swapped

             data_out(i  ,j  ,ID) = data_out(i  ,j  ,ID) + flux_y(ID)*dtdx
             data_out(i  ,j  ,IP) = data_out(i  ,j  ,IP) + flux_y(IP)*dtdx
             data_out(i  ,j  ,IU) = data_out(i  ,j  ,IU) + flux_y(IV)*dtdx ! watchout IU and IV swapped
             data_out(i  ,j  ,IV) = data_out(i  ,j  ,IV) + flux_y(IU)*dtdx ! watchout IU and IV swapped

          end do ! end do j
       end do ! end do i

    else if (implementationVersion .eq. 1) then

       ! convert conservative variable into primitives ones for the entire domain
       call convertToPrimitives(data_in)

       ! trace computation: fill global arrays qm_x, qm_y, qp_x, qp_y
       call computeTrace(dt)

       ! Compute flux via Riemann solver and update (time integration)
       call computeFluxesAndUpdate(data_out,dt)

    end if

    call timerStop(godunov_timer)

  end subroutine godunov_unsplit_cpu

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert conservative variables array into primitive var array (q)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine convertToPrimitives(data)
    implicit none

    ! dummy variables
    real(fp_kind), dimension(isize, jsize, nbVar), intent(inout) :: data

    ! local variables
    ! primitive variable state vector
    real(fp_kind), dimension(nbVar) :: qLoc
    real(fp_kind)                   :: c
    integer :: i,j

    do j=1,jsize
       do i=1,isize

          call computePrimitives(data, i, j, c, qLoc)

          ! copy q state in q global
          q(i,j,:) = qLoc

       end do ! end do i
    end do ! end do j

  end subroutine convertToPrimitives

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute trace (only used in implementation version 1), i.e.
  !! fill global array qm_x, qmy, qp_x, qp_y
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine computeTrace(dt)
    implicit none

    ! dummy variables
    real(fp_kind), intent(in)    :: dt

    ! local variables
    integer :: i,j,iVar
    real(fp_kind), dimension(nbVar)   :: qLoc ! local primitive variables
    real(fp_kind), dimension(nbVar)   :: qPlusX
    real(fp_kind), dimension(nbVar)   :: qMinusX
    real(fp_kind), dimension(nbVar)   :: qPlusY
    real(fp_kind), dimension(nbVar)   :: qMinusY
    real(fp_kind), dimension(2,nbVar) :: dq
    real(fp_kind), dimension(2,nbVar) :: qm, qp
    real(fp_kind) :: dtdx
    real(fp_kind) :: dtdy

    dtdx = dt / dx
    dtdy = dt / dy

    do j=2,jsize-1
       do i=2,isize-1

          ! get primitive variables state vector
          do iVar=1,nbVar
             qLoc   (iVar) = q(i  ,j  ,iVar)
             qPlusX (iVar) = q(i+1,j  ,iVar)
             qMinusX(iVar) = q(i-1,j  ,iVar)
             qPlusY (iVar) = q(i  ,j+1,iVar)
             qMinusY(iVar) = q(i  ,j-1,iVar)
          end do

	  ! get hydro slopes dq
          call slope_unsplit_hydro_2d(qLoc, qPlusX, qMinusX, qPlusY, qMinusY, dq)

          ! compute qm, qp
          call trace_unsplit_hydro_2d(qLoc, dq, dtdx, dtdy, qm, qp)

          ! store qm, qp : only what is really needed
          do iVar=1,nbVar
             qm_x(i,j,iVar) = qm(1,iVar)
             qp_x(i,j,iVar) = qp(1,iVar)
             qm_y(i,j,iVar) = qm(2,iVar)
             qp_y(i,j,iVar) = qp(2,iVar)
          end do ! end do iVar

       end do ! end do i
    end do ! end do j

  end subroutine computeTrace

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute flux via Riemann solver and update (time integration)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine computeFluxesAndUpdate(data,dt)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(isize, jsize, nbVar), intent(inout) :: data
    real(fp_kind), intent(in)    :: dt

    ! local variables
    integer :: i,j
    real(fp_kind), dimension(nbVar)     :: qleft, qright
    real(fp_kind), dimension(nbVar)     :: flux_x, flux_y
    real(fp_kind), dimension(nbVar)     :: qgdnv

    real(fp_kind) :: dtdx
    real(fp_kind) :: dtdy

    dtdx = dt / dx
    dtdy = dt / dy
    
    flux_x = 0
    flux_y = 0

    do j=ghostWidth+1,jsize-ghostWidth+1
       do i=ghostWidth+1,isize-ghostWidth+1
          
          !
          ! Solve Riemann problem at X-interfaces and compute
          ! X-fluxes
          !
          qleft(ID)   = qm_x(i-1,j,ID)
          qleft(IP)   = qm_x(i-1,j,IP)
          qleft(IU)   = qm_x(i-1,j,IU)
          qleft(IV)   = qm_x(i-1,j,IV)

          qright(ID)  = qp_x(i  ,j,ID)
          qright(IP)  = qp_x(i  ,j,IP)
          qright(IU)  = qp_x(i  ,j,IU)
          qright(IV)  = qp_x(i  ,j,IV)

          ! compute hydro flux_x
          call riemann_2d(qleft,qright,qgdnv,flux_x)

	  !
          ! Solve Riemann problem at Y-interfaces and compute Y-fluxes
          !
          qleft(ID)   = qm_y(i,j-1,ID)
          qleft(IP)   = qm_y(i,j-1,IP)
          qleft(IU)   = qm_y(i,j-1,IV) ! watchout IU, IV permutation
          qleft(IV)   = qm_y(i,j-1,IU) ! watchout IU, IV permutation

          qright(ID)  = qp_y(i,j  ,ID)
          qright(IP)  = qp_y(i,j  ,IP)
          qright(IU)  = qp_y(i,j  ,IV) ! watchout IU, IV permutation
          qright(IV)  = qp_y(i,j  ,IU) ! watchout IU, IV permutation

          ! compute hydro flux_y
          call riemann_2d(qleft,qright,qgdnv,flux_y)

          !
          ! update hydro array
          !
          data(i-1,j  ,ID) = data(i-1,j  ,ID) - flux_x(ID)*dtdx
          data(i-1,j  ,IP) = data(i-1,j  ,IP) - flux_x(IP)*dtdx
          data(i-1,j  ,IU) = data(i-1,j  ,IU) - flux_x(IU)*dtdx
          data(i-1,j  ,IV) = data(i-1,j  ,IV) - flux_x(IV)*dtdx

          data(i  ,j  ,ID) = data(i  ,j  ,ID) + flux_x(ID)*dtdx
          data(i  ,j  ,IP) = data(i  ,j  ,IP) + flux_x(IP)*dtdx
          data(i  ,j  ,IU) = data(i  ,j  ,IU) + flux_x(IU)*dtdx
          data(i  ,j  ,IV) = data(i  ,j  ,IV) + flux_x(IV)*dtdx

          data(i  ,j-1,ID) = data(i  ,j-1,ID) - flux_y(ID)*dtdy
          data(i  ,j-1,IP) = data(i  ,j-1,IP) - flux_y(IP)*dtdy
          data(i  ,j-1,IU) = data(i  ,j-1,IU) - flux_y(IV)*dtdy ! watchout IU and IV swapped
          data(i  ,j-1,IV) = data(i  ,j-1,IV) - flux_y(IU)*dtdy ! watchout IU and IV swapped

          data(i  ,j  ,ID) = data(i  ,j  ,ID) + flux_y(ID)*dtdy
          data(i  ,j  ,IP) = data(i  ,j  ,IP) + flux_y(IP)*dtdy
          data(i  ,j  ,IU) = data(i  ,j  ,IU) + flux_y(IV)*dtdy ! watchout IU and IV swapped
          data(i  ,j  ,IV) = data(i  ,j  ,IV) + flux_y(IU)*dtdy ! watchout IU and IV swapped

       end do ! end do i
    end do ! end do j

  end subroutine computeFluxesAndUpdate

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Hydrodynamical Implosion Test :
  !! http://www.astro.princeton.edu/~jstone/tests/implode/Implode.html
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_implode(data)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(isize, jsize, nbVar), intent(inout) :: data

    ! local variables
    integer :: i,j
    real(fp_kind) :: tmp

    do j=1,jsize
       do i=1,isize
          tmp = 1.0*(i+decx-ghostWidth-1)/nx + 1.0*(j+decy-ghostWidth-1)/ny
          if (tmp .gt. 0.5) then
            
             data(i,j,ID)=1.0+tmp
             data(i,j,IP)=1.0/(gamma0-1.0)
             data(i,j,IU)=0.0
             data(i,j,IV)=0.0
          else
             data(i,j,ID)=0.125+tmp
             data(i,j,IP)=0.14/(gamma0-1.0)
             data(i,j,IU)=0.0
             data(i,j,IV)=0.0               
          end if
       end do
    end do

  end subroutine init_implode

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! output routine (VTK file format, ASCII, VtkImageData)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine saveVTK(data,iStep)
    implicit none
    ! dummy variables
    real   (fp_kind), dimension(isize_tot, jsize_tot, nbVar), intent(inout) :: data
    !real   (fp_kind), dimension(isize, jsize, nbVar), intent(inout) :: data
    integer(int_kind) :: iStep
    !integer, intent(in) :: myRank

    ! local variables
    integer :: j,iVar
    integer :: error
    character(LEN=80) :: filename
    character(LEN=8)  :: filenum
    character(1), parameter :: endl  = char(10)  ! end of line
    character(500) :: charBuf
    character(7) :: floatType

    !if(myRank==0) then

    if (useDoublePrecision()) then
       write(floatType,'(a)') 'Float64'
    else
       write(floatType,'(a)') 'Float32'
    end if
    write (filenum,'(i8.8)') iStep
    filename='euler2d_' // filenum // '.vti'

    !open(10,file=filename,status='replace',form='unformatted',action='write',iostat=error)
    !open(10,file=filename,status='replace',access='stream',action='write',iostat=error)
    open(10,file=filename,iostat=error)

    ! write header
    write(10,'(a)') '<?xml version="1.0"?>'//endl
    if (isBigEndian()) then
       write(10,'(a)') '<VTKFile type="ImageData" version="0.1" byte_order="BigEndian">'//endl
    else
       write(10,'(a)') '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//endl
    end if

    ! write mesh extent
    write(charBuf,fmt='(6(I7))',iostat=error) 1,isize_tot-2*ghostWidth+1,1,jsize_tot-2*ghostWidth+1,1,2
    !write(charBuf,fmt='(6(I7))',iostat=error) 1,isize-2*ghostWidth+1,1,jsize-2*ghostWidth+1,1,2
    write(10,'(a)') repeat(' ',2)//'<ImageData WholeExtent="'//trim(charBuf)//'"'
    write(10,'(a)') ' Origin="0 0 0" Spacing="1 1 1">'//endl
    write(10,'(a)') repeat(' ',2)//'<Piece Extent="'//trim(charBuf)//'">'//endl

    write(10,'(a)') repeat(' ',3)//'<PointData>'//endl
    write(10,'(a)') repeat(' ',4)//'</PointData>'//endl
    write(10,'(a)') repeat(' ',4)//'<CellData>'//endl

    ! write data array (ascii), remove ghost cells
    do iVar=1,nbVar
       write(10,'(a)') repeat(' ',4)//'<DataArray type="'//trim(floatType)// &
            & '" Name="'//varNames(iVar)//'" format="ascii" >'//endl

       do j=ghostWidth+1,jsize_tot-ghostWidth
          write(10,*) data(ghostWidth+1:isize_tot-ghostWidth,j,iVar)
!       do j=ghostWidth+1,jsize-ghostWidth
!          write(10,*) data(ghostWidth+1:isize-ghostWidth,j,iVar)
       end do

       write(10,'(a)') repeat(' ',4)//'</DataArray>'//endl
    end do

    ! write footer
    write(10,'(a)') repeat(' ',4)//'</CellData>'//endl
    write(10,'(a)') repeat(' ',2)//'</Piece>'//endl
    write(10,'(a)') repeat(' ',2)//'</ImageData>'//endl
    write(10,'(a)') '</VTKFile>'//endl

    close(10)
 
    !end if

  end subroutine saveVTK

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Fill ghost cells according to border condition :
  !! absorbant, reflexive or periodic
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_boundaries(data)
    implicit none
    ! dummy variables
    real(fp_kind), dimension(isize, jsize, nbVar), intent(inout) :: data

    ! local variables
    integer ::i,j,i0,j0,iVar
    real(fp_kind) :: sign

    ! write inner border conditions in communication arrays
    call updateS(data)
    ! send and receive communication arrays
    call comm

    ! boundary xmin
    if(coord_x==0) then
     do iVar=1,nbVar
        do i=1,ghostWidth
           sign=1.0
           if(boundary_type_xmin==1)then
              i0=2*ghostWidth+1-i
              if(iVar==IU)sign=-1.0
           else if(boundary_type_xmin==2)then
              i0=ghostWidth+1
           else ! periodic
             i0=nx+i
          end if
          do j=ghostWidth+1,jsize-ghostWidth
             data(i,j,iVar)=data(i0,j,iVar)*sign
          end do
       end do
    end do
    else ! inner cell
      do iVar=1,nbVar
        do i=1,ghostWidth
          do j=1+ghostWidth,jsize-ghostWidth
             data(i,j,iVar)=leftR(i,j-ghostWidth,iVar)
          end do
        end do
      end do
    end if

    ! boundary xmax
    if(coord_x==size_x_max-1) then
    do iVar=1,nbVar
       do i=1+isize-ghostWidth,isize
          sign=1.0
          if(boundary_type_xmax==1)then
             i0=2*isize-2*ghostWidth+1-i
             if(iVar==IU)sign=-1.0
          else if(boundary_type_xmax==2)then
             i0=isize-ghostWidth
          else ! periodic
             i0=i-nx
          end if
          do j=ghostWidth+1,jsize-ghostWidth
             data(i,j,iVar)=data(i0,j,iVar)*sign
          end do
       end do
    end do
    else ! inner cell
       do iVar=1,nbVar
          do i=1+isize-ghostWidth,isize
             do j=ghostWidth+1,jsize-ghostWidth
                data(i,j,iVar)=rightR(i-isize+ghostWidth,j-ghostWidth,iVar)
             end do
          end do
       end do
    end if

    ! boundary ymin
    if(coord_y==0) then
    do iVar=1,nbVar
       do j=1,ghostWidth
          sign=1.0
          if(boundary_type_ymin==1)then
             j0=2*ghostWidth+1-j
             if(iVar==IV)sign=-1.0
          else if(boundary_type_ymin==2)then
             j0=ghostWidth+1
          else ! periodic
             j0=ny+j
          end if
          do i=ghostWidth+1,isize-ghostWidth
             data(i,j,iVar)=data(i,j0,iVar)*sign
          end do
       end do
    end do
    else ! inner cell
       do iVar=1,nbVar
          do j=1,ghostWidth
             do i=1+ghostwidth,isize-ghostWidth
                data(i,j,iVar)=bottomR(i-ghostWidth,j,iVar)
             end do
          end do
       end do
    end if

    ! boundary ymax
    if (coord_y == size_y_max-1) then
    do iVar=1,nbVar
       do j=jsize-ghostWidth+1,jsize
          sign=1.0
          if(boundary_type_ymax==1)then
             j0=2*jsize-2*ghostWidth+1-j
             if(iVar==IV)sign=-1.0
          else if(boundary_type_ymax==2)then
             j0=ny+ghostWidth
          else ! periodic
             j0=j-ny
          end if
          do i=ghostWidth+1,isize-ghostWidth
             data(i,j,iVar)=data(i,j0,iVar)*sign
          end do
       end do
    end do
    else ! inner cell
     do iVar=1,nbVar
      do j = jsize-ghostWidth+1,jsize
        do i = ghostWidth+1,isize-ghostWidth
            data(i,j,iVar)=topR(i-ghostWidth,j-jsize+ghostWidth,iVar)
        end do
      end do
    end do
   end if 

  end subroutine make_boundaries

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! create the global array for printing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine reconstitute(myRank, nbTask)

    implicit none
    
    ! dummy variables
    integer, intent(in) :: myRank, nbTask
    ! local vairables
    integer :: ierr, i, j, iVar
    integer, dimension(MPI_STATUS_SIZE) :: status

    ! processor 0 sends the global_array
    if(myRank == 0) then
      call MPI_SEND(u_tot, isize_tot*jsize_tot*nbVar, MPI_REAL, 1, 1, MPI_COMM_WORLD, ierr)
    end if
 
    ! each processor writes its data in it
    if(myRank > 0) then
      call MPI_RECV(u_tot, isize_tot*jsize_tot*nbVar, MPI_REAL, myRank-1, 1, MPI_COMM_WORLD, status, ierr)
      do iVar = 1,nbVar
        do i = 1+ghostWidth,isize-ghostWidth
          do j = 1+ghostWidth,jsize-ghostWidth
            u_tot(i+decx,j+decy,iVar) = u(i,j,iVar)
          end do
        end do
      end do
      call MPI_SEND(u_tot, isize_tot*jsize_tot*nbVar, MPI_REAL, modulo(myRank+1,nbTask), 1, MPI_COMM_WORLD, ierr) 
    end if

    ! processor 0 receives the completed global array and writes in it
    if(myRank == 0) then
      call MPI_RECV(u_tot, isize_tot*jsize_tot*nbVar, MPI_REAL, nbTask-1, 1, MPI_COMM_WORLD, status, ierr)
      do iVar = 1,nbVar
        do i = 1+ghostWidth,isize-ghostWidth
          do j = 1+ghostWidth,jsize-ghostWidth
            u_tot(i+decx,j+decy,iVar) = u(i,j,iVar)
          end do
        end do
      end do
    end if

    end subroutine

end module HydroRun
