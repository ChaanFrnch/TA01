!!!! -*- Mode: F90 -*- !!!!
!> \file main.f90
!> \brief 2D Euler solver

program euler2d

  use HydroParameters  ! get routines initHydroParameters, printHydroParameters
  use HydroRun         ! get computing routines and utilities (init, boundaries, ...)
  use Monitoring       ! get timer routines
  use Partitioner      ! partition of the domain
  use Communication    ! for parallel communications
  use mpi              ! for parallelization

  implicit none

  real   (fp_kind)  :: t=0
  real   (fp_kind)  :: dt=0
  integer :: nbTask, myRank, ierr

  call MPI_Init(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr) ! get number of processors
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr) ! get the rank of each processor

  call initHydroParameters(nbTask, myRank)
  if (myRank == 0 ) then
    call printHydroParameters()
  end if

  ! init domain
  call initHydroRun
  call compute_dt( dt, modulo(nStep,2) ) ! most constraining dt from all processors 

  ! init boundaries
  call initS
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! start computation
  call timerStart(total_timer)

  ! main loop
  do while (t < tEnd .and. nStep < nStepmax) ! loop over time and maximum number of steps
     ! output
     if ( modulo(nStep,nOutput) == 0) then ! print every nOutput
        call reconstitute(myRank, nbTask)  ! assemble the total array
        if (myRank == 0) then
        write(*,*) 'Output results at step ',nStep, 'dt ',dt
        call timerStart(io_timer)
        call saveVTK(u_tot,nStep)
        call timerStop(io_timer)
        end if
     end if

     ! compute dt_min
     call compute_dt( dt, modulo(nStep,2) ) ! most constraining dt from all processors 

     ! perform one step integration
     call godunov_unsplit(dt)
     nStep = nStep+1
  end do

  ! end of computation
  call timerStop(total_timer)
  call cleanupHydroRun()

  ! print monitoring
 if (myRank == 0) then
  write(*,*) 'total      time : ', total_timer%elapsed,     'secondes'
  write(*,*) 'compute    time : ', godunov_timer%elapsed,   'secondes', 100*godunov_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'io         time : ', io_timer%elapsed,        'secondes', 100*io_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'boundaries time : ', boundaries_timer%elapsed,'secondes', 100*boundaries_timer%elapsed / total_timer%elapsed, '%'

  write(*,*) 'Perf             : ',nStep*isize*jsize/total_timer%elapsed, ' number of cell-updates/s'
 end if

  call MPI_Finalize(ierr)

end program euler2d
