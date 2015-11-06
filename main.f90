!!!! -*- Mode: F90 -*- !!!!
!> \file main.f90
!> \brief 2D Euler solver

program euler2d

  use HydroParameters  ! get routines initHydroParameters, printHydroParameters
  use HydroRun         ! get computing routines and utilities (init, boundaries, ...)
  use Monitoring       ! get timer routines
  use Partitioner      ! partition of the domain
  use mpi              ! for parallelization

  implicit none

  real   (fp_kind)  :: t=0
  real   (fp_kind)  :: dt=0
  !integer :: nbpowers, nbcores
  integer :: nbTask, myRank, ierr
  !integer,dimension(:),allocatable :: sizes_x, sizes_y

  call MPI_Init(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

  call initHydroParameters(nbTask, myRank, sizes_x, sizes_y)
  call printHydroParameters()

  ! init domain
  call initHydroRun()
  call compute_dt( dt, modulo(nStep,2) ) ! pas adaptatif, calcule a chaque fois pour stabilite
  write(*,*) 'Initial value for dt ',dt

  ! init boundaries
  call make_boundaries(u)

  ! start computation
  write(*,*) 'Start computation....'
  call timerStart(total_timer)

  ! main loop
  do while (t < tEnd .and. nStep < nStepmax) ! boucle sur le temps et le nb de pas
     ! output
     if ( modulo(nStep,nOutput) == 0) then ! impression tous les nOutput
        write(*,*) 'Output results at step ',nStep, 'dt ',dt
        call timerStart(io_timer)
        call saveVTK(u,nStep)
        call timerStop(io_timer)
     end if
     
     ! compute dt
     call compute_dt( dt, modulo(nStep,2) )

     ! perform one step integration
     call godunov_unsplit(dt)

     nStep = nStep+1

  end do

  ! end of computation
  call timerStop(total_timer)

  call cleanupHydroRun()

  ! print monitoring
  write(*,*) 'total      time : ', total_timer%elapsed,     'secondes'
  write(*,*) 'compute    time : ', godunov_timer%elapsed,   'secondes', 100*godunov_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'io         time : ', io_timer%elapsed,        'secondes', 100*io_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'boundaries time : ', boundaries_timer%elapsed,'secondes', 100*boundaries_timer%elapsed / total_timer%elapsed, '%'

  write(*,*) 'Perf             : ',nStep*isize*jsize/total_timer%elapsed, ' number of cell-updates/s'

  call MPI_Finalize(ierr)

end program euler2d
