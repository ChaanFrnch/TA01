module Communication

   use HydroParameters
   use Monitoring
   use mpi

   real(fp_kind), dimension(:,:,:), allocatable :: topS, bottomS, leftS, rightS, topR, bottomR, leftR, rightR

contains

   subroutine initS(data)

   implicit none

   integer :: i, j, iVar
   real(fp_kind), dimension(isize, jsize, nbVar), intent(in) :: data

   allocate (topS(ghostWidth, jsize-2*ghostWidth, nbVar))
   allocate (bottomS(ghostWidth, jsize-2*ghostWidth, nbVar))
   allocate (leftS(isize-2*ghostWidth, ghostWidth, nbVar)) 
   allocate (rightS(isize-2*ghostWidth, ghostWidth, nbVar))
   allocate (topR(ghostWidth, jsize-2*ghostWidth, nbVar))
   allocate (bottomR(ghostWidth, jsize-2*ghostWidth, nbVar))
   allocate (leftR(isize-2*ghostWidth, ghostWidth, nbVar)) 
   allocate (rightR(isize-2*ghostWidth, ghostWidth, nbVar))

       do iVar = 1, nbVar
       do i = 1, ghostWidth
          do j = 1, jsize-2*ghostWidth
              topS(i,j,iVar) = data(i+ghostWidth,j+ghostWidth,iVar)
              bottomS(i,j,iVar) = data(isize-ghostWidth+i, j+ghostWidth, iVar)
          end do
       end do

       do j = 1, ghostWidth
          do i = 1, isize-2*ghostWidth
              leftS(i,j,iVar) = data(i+ghostWidth, j+ghostWidth, iVar)
              rightS(i,j,iVar) = data(i+ghostWidth, j+jsize-ghostWidth, iVar)
          end do
       end do
       end do

   end subroutine initS

   subroutine comm

      if(coord_y>1) then
         call MPI_ISEND(topS, 1, MPI_REAL, coord_y-1+size_x_max*coord_x, 1, MPI_COMM_WORLD, ierr)
         call MPI_IRECV(topR, 1, MPI_REAL, coord_y-1+size_x_max*coord_x, 2, MPI_COMM_WORLD, ierr)
      end if

      if(coord_y<size_y_max) then
         call MPI_ISEND(bottomS, 1, MPI_REAL, coord_y+1+size_x_max*coord_x, 2, MPI_COMM_WORLD, ierr)
         call MPI_IRECV(bottomR, 1, MPI_REAL, coord_y+1+size_x_max*coord_x, 1, MPI_COMM_WORLD, ierr)
      end if

      if(coord_x>1) then
         call MPI_ISEND(leftS, 1, MPI_REAL, coord_y+size_x_max*(coord_x-1), 3, MPI_COMM_WORLD, ierr)
         call MPI_IRECV(leftR, 1, MPI_REAL, coord_y+size_x_max*(coord_x-1), 4, MPI_COMM_WORLD, ierr)
      end if

      if(coord_x<size_x_max) then
         call MPI_ISEND(rightS, 1, MPI_REAL, coord_y+size_x_max*(coord_x+1), 4, MPI_COMM_WORLD, ierr)
         call MPI_IRECV(rightR, 1, MPI_REAL, coord_y+size_x_max*(coord_x+1), 3, MPI_COMM_WORLD, ierr)
      end if

   end subroutine comm

end module Communication
