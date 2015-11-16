module Communication

   use HydroParameters
   use Monitoring
   use mpi

   real(fp_kind), dimension(:,:,:), allocatable :: topS, bottomS, leftS, rightS, topR, bottomR, leftR, rightR
   integer :: size_top,size_right

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

   integer :: i, j, iVar
   real(fp_kind), dimension(isize, jsize, nbVar), intent(in) :: data

   do iVar = 1, nbVar
       do i = 1,ghostWidth
          do j = 1,jsize-2*ghostWidth
              topS(i,j,iVar) = data(i+ghostWidth, j+ghostWidth,iVar)
              bottomS(i,j,iVar) = data(i+isize-2*ghostWidth, j+ghostWidth, iVar)
          end do
       end do

       do i = 1, isize-2*ghostWidth
          do j = 1,ghostWidth
              leftS(i,j,iVar) = data(i+ghostWidth,j+ghostWidth , iVar)
              rightS(i,j,iVar) = data(i+ghostWidth, j+jsize-2*ghostWidth, iVar)
          end do
         end do
   end do

  size_top = (jsize-2*ghostwidth)*ghostWidth*nbVar
  size_right = (isize-2*ghostwidth)*ghostWidth*nbVar

  end subroutine updateS

  subroutine comm

  integer :: nb_loop, index_loop,index_x,index_y

!write(*,*) size_x_max , " " , size_y_max

  nb_loop = size_x_max

!top send
  do index_loop = 0,nb_loop
  index_x = size_x_max - index_loop
   if (coord_x < size_x_max-1 .AND. coord_x == index_x) then
    call MPI_RECV(bottomR, size_top, MPI_REAL, coord_y+size_x_max*(index_x+1), 1, MPI_COMM_WORLD, mpistat, ierr)
   end if 
   if (coord_x > 0 .AND. coord_x == index_x ) then
    call MPI_SEND(topS, size_top, MPI_REAL, coord_y + size_x_max*(index_x-1), 1, MPI_COMM_WORLD, mpistat, ierr) 
   end if
  end do

!bot send
  do index_loop = 0,nb_loop
   if (coord_x > 0 .AND. coord_x == index_loop) then
    call MPI_RECV(topR, size_top, MPI_REAL, coord_y+size_x_max*(index_loop-1), 1, MPI_COMM_WORLD, mpistat, ierr)
   end if 
   if (coord_x < size_x_max-1 .AND. coord_x == index_loop ) then
    call MPI_SEND(bottomS, size_top, MPI_REAL, coord_y + size_x_max*(index_loop+1), 1, MPI_COMM_WORLD, mpistat, ierr) 
   end if
  end do

  nb_loop = size_y_max 


  !left send
  do index_loop = 0,nb_loop
    index_y = size_y_max - index_loop
   if (coord_y < size_y_max-1 .AND. coord_y == index_y) then
      call MPI_RECV(rightR, size_right, MPI_REAL, size_x_max*coord_x + (index_y +1), 1, MPI_COMM_WORLD, mpistat, ierr)
    end if 
    if (coord_y > 0 .AND. coord_y == index_y ) then
      call MPI_SEND(leftS, size_right, MPI_REAL, (index_y-1) + size_x_max*coord_x, 1, MPI_COMM_WORLD, mpistat, ierr) 
    end if
  end do

!right send
  do index_loop = 0,nb_loop
    if (coord_y > 0 .AND. coord_y == index_loop) then
     call MPI_RECV(leftR, size_right, MPI_REAL, (index_loop-1)+size_x_max*coord_x, 1, MPI_COMM_WORLD, mpistat, ierr)
    end if 
    if (coord_y < size_y_max-1 .AND. coord_y == index_loop ) then
     call MPI_SEND(rightS, size_right, MPI_REAL, (index_loop+1)+size_x_max*coord_x, 1, MPI_COMM_WORLD, mpistat, ierr) 
    end if
  end do


   end subroutine comm

end module Communication
