module Monitoring

  use HydroPrecision

  ! declare a Timer type
  type Timer
     real(fp_kind) :: elapsed=0.0
     real(fp_kind) :: start=0.0
     real(fp_kind) :: stop=0.0
  end type Timer
  
  type(Timer) :: total_timer
  type(Timer) :: godunov_timer   
  type(Timer) :: boundaries_timer
  type(Timer) :: io_timer

  contains

    ! start timer
    subroutine timerStart(t)
      implicit none
      type(Timer), intent(inout) :: t

      call cpu_time(t%start)

    end subroutine timerStart

    ! stop timer and accumulate timings
    subroutine timerStop(t)
      implicit none
      type(Timer), intent(inout) :: t

      call cpu_time(t%stop)
      t%elapsed = t%elapsed + t%stop - t%start

    end subroutine timerStop

end module Monitoring
