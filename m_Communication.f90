module Communication

   use HydroParameters
   use Monitoring
   use mpi

   real(fp_kind), dimension(:,:,:), allocatable :: topS, bottomS, leftS, rightS, topR, bottomR, leftR, rightR

contains

   subroutine initS

   allocate (topS(ghostWidth,jsize-2*ghostWidth,nbVar))
   allocate (topR(ghostWidth,jsize-2*ghostWidth,nbVar))

   allocate (bottomS(ghostWidth,jsize-2*ghostWidth,nbVar))
   allocate (bottomR(ghostWidth,jsize-2*ghostWidth,nbVar))

   allocate (leftS(isize-2*ghostWidth,ghostWidth,nbVar))
   allocate (leftR(isize-2*ghostWidth,ghostWidth,nbVar)) 
 
   allocate (rightS(isize-2*ghostWidth,ghostWidth,nbVar))
   allocate (rightR(isize-2*ghostWidth,ghostWidth,nbVar))

   end subroutine initS

   subroutine updateS(data)

   implicit none 

   integer :: i, j, iVar,size_top,size_right
   real(fp_kind), dimension(isize, jsize, nbVar), intent(in) :: data



   do iVar = 1, nbVar
       do i = 1,ghostWidth
          do j = 1,jsize-2*ghostWidth
              topS(i,j,iVar) = data(i+ghostWidth, j+ghostWidth,iVar)
              bottomS(i,j,iVar) = data(i+ghostWidth, j+jsize-2*ghostWidth, iVar)
          end do
       end do

       do i = 1, jsize-2*ghostWidth
          do j = 1,ghostWidth
              leftS(i,j,iVar) = data(i+ghostWidth, j+ghostWidth, iVar)
              rightS(i,j,iVar) = data(i+isize-2*ghostWidth, j+ghostWidth, iVar)
          end do
         end do
   end do

  size_top = (jsize-2*ghostwidth)*ghostWidth*nbVar
  !write(*,*) 'size_top = ', size_top
  !write(*,*) 'la taille du tableau est ', size(topS)
  size_right = (isize-2*ghostwidth)*ghostWidth*nbVar

  end subroutine updateS

  subroutine comm

      if(coord_y>0) then
         write(*,*) 'I am proc ', coord_x, coord_y, 'ok1'
         call MPI_ISEND(topS, size_top, MPI_REAL, coord_y-1+size_x_max*coord_x, 1, MPI_COMM_WORLD, ierr)
         write(*,*) 'ok2'
         call MPI_IRECV(topR, size_top, MPI_REAL, coord_y-1+size_x_max*coord_x, 2, MPI_COMM_WORLD, ierr) 
      end if

      if(coord_y<size_y_max-1) then
         write(*,*) 'ok3'
         call MPI_ISEND(bottomS, size_top, MPI_REAL, coord_y+1+size_x_max*coord_x, 2, MPI_COMM_WORLD, ierr)
         write(*,*) 'ok4'
         call MPI_IRECV(bottomR, size_top, MPI_REAL, coord_y+1+size_x_max*coord_x, 1, MPI_COMM_WORLD, ierr)
      end if

      if(coord_x>0) then
         write(*,*) 'ok5'
         call MPI_ISEND(leftS, size_right, MPI_REAL, coord_y+size_x_max*(coord_x-1), 3, MPI_COMM_WORLD, ierr)
         write(*,*) 'ok6'
         call MPI_IRECV(leftR, size_right, MPI_REAL, coord_y+size_x_max*(coord_x-1), 4, MPI_COMM_WORLD, ierr)
      end if

      if(coord_x<size_x_max-1) then
         write(*,*) 'ok7'
         call MPI_ISEND(rightS, size_right, MPI_REAL, coord_y+size_x_max*(coord_x+1), 4, MPI_COMM_WORLD, ierr)
         write(*,*) 'ok8'
         call MPI_IRECV(rightR, size_right, MPI_REAL, coord_y+size_x_max*(coord_x+1), 3, MPI_COMM_WORLD, ierr)
      end if

   end subroutine comm

end module Communication
