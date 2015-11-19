!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!< Defines all variables needed to perform parallel communications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Communication

   use HydroParameters
   use Monitoring
   use mpi

   ! communication parameters
   real(fp_kind), dimension(:,:,:), allocatable :: topS, bottomS, leftS, rightS, topR, bottomR, leftR, rightR
   integer :: size_top,size_right

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! initialization of the communication arrays
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine initS

   ! top array (to send and to receive)
   allocate (topS(isize-2*ghostWidth,ghostWidth,nbVar))
   allocate (topR(isize-2*ghostWidth,ghostWidth,nbVar))

   ! bottom array (to send and to receive)
   allocate (bottomS(isize-2*ghostWidth,ghostWidth,nbVar))
   allocate (bottomR(isize-2*ghostWidth,ghostWidth,nbVar))

   ! left array (to send and to receive)
   allocate (leftS(ghostWidth,jsize-2*ghostWidth,nbVar))
   allocate (leftR(ghostWidth,jsize-2*ghostWidth,nbVar)) 
 
   ! right array (to send and to receive)
   allocate (rightS(ghostWidth,jsize-2*ghostWidth,nbVar))
   allocate (rightR(ghostWidth,jsize-2*ghostWidth,nbVar))

   end subroutine initS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! write the boundaries values in the communication arrays
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine updateS(data)

   implicit none 

   ! local variables
   integer :: i, j, iVar
   real(fp_kind), dimension(isize, jsize, nbVar), intent(in) :: data

   do iVar = 1, nbVar

       do i = 1,isize-2*ghostWidth
         do j = 1,ghostWidth
           topS(i,j,iVar)   = data(i+ghostWidth, j+jsize-2*ghostWidth,iVar)
           bottomS(i,j,iVar)= data(i+ghostWidth, j+ghostWidth, iVar)
         end do
       end do

       do j = 1,jsize-2*ghostWidth
         do i = 1,ghostWidth
           rightS(i,j,iVar) = data(i+isize-2*ghostWidth, j+ghostWidth, iVar)
           leftS(i,j,iVar)  = data(i+ghostWidth,j+ghostWidth , iVar)
         end do
       end do
   end do

  ! sizes of the communication arrays
  size_top   = (isize-2*ghostwidth)*ghostWidth*nbVar
  size_right = (jsize-2*ghostwidth)*ghostWidth*nbVar

  end subroutine updateS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! sending and receiving communication arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine comm

  ! local variables
  integer :: nb_loop, index_loop,index_x,index_y

  nb_loop = size_y_max-1

  !top send
  do index_loop = 0,nb_loop
   index_y = index_loop
   index_x = coord_x
   if (coord_y > 0 .AND. coord_y == index_y) then
    call MPI_RECV(bottomR, size_top, MPI_REAL, index_x+size_x_max*(index_y-1) ,1 , MPI_COMM_WORLD, mpistat, ierr)
   end if 
   if (coord_y < nb_loop .AND. coord_y == index_y ) then
    call MPI_SEND(topS   , size_top, MPI_REAL, index_x+size_x_max*(index_y+1) ,1 , MPI_COMM_WORLD, mpistat, ierr)
   end if
  end do

  !bottom send
  do index_loop = 0,nb_loop
  index_y = nb_loop - index_loop
  index_x = coord_x
   if (coord_y < nb_loop .AND. coord_y == index_y) then
    call MPI_RECV(topR   , size_top, MPI_REAL, index_x+size_x_max*(index_y+1), 2, MPI_COMM_WORLD, mpistat, ierr)
   end if 
   if (coord_y > 0 .AND. coord_y == index_y ) then
    call MPI_SEND(bottomS, size_top, MPI_REAL, index_x+size_x_max*(index_y-1), 2, MPI_COMM_WORLD, mpistat, ierr) 
   end if
  end do

  nb_loop = size_x_max-1 

  !right send
  do index_loop = 0,nb_loop
   index_x = index_loop
   index_y = coord_y
   if (coord_x > 0 .AND. coord_x == index_x) then
    call MPI_RECV(leftR , size_top, MPI_REAL, (index_x-1)+size_x_max*index_y , 3, MPI_COMM_WORLD, mpistat, ierr)
   end if 
   if (coord_x < nb_loop .AND. coord_x == index_x ) then
    call MPI_SEND(rightS, size_top, MPI_REAL, (index_x+1)+size_x_max*index_y , 3, MPI_COMM_WORLD, mpistat, ierr)
   end if
  end do

  !bottom send
  do index_loop = 0,nb_loop
  index_x = nb_loop - index_loop
  index_y = coord_y
   if (coord_x < nb_loop .AND. coord_x == index_x) then
    call MPI_RECV(rightR, size_top, MPI_REAL, (index_x+1)+size_x_max*index_y, 4, MPI_COMM_WORLD, mpistat, ierr)
   end if 
   if (coord_x > 0 .AND. coord_x == index_x ) then
    call MPI_SEND(leftS , size_top, MPI_REAL, (index_x-1)+size_x_max*index_y, 4, MPI_COMM_WORLD, mpistat, ierr) 
   end if
  end do

  end subroutine comm

end module Communication
