module Communication

   use HydroParameters
   use Monitoring
   use mpi

   real(fp_kind), dimension(:,:,:), allocatable :: topS, bottomS, leftS, rightS, topR, bottomR, leftR, rightR
   integer :: size_top,size_right

contains

   subroutine initS

   allocate (topS(isize-2*ghostWidth,ghostWidth,nbVar))
   allocate (topR(isize-2*ghostWidth,ghostWidth,nbVar))

   allocate (bottomS(isize-2*ghostWidth,ghostWidth,nbVar))
   allocate (bottomR(isize-2*ghostWidth,ghostWidth,nbVar))

   allocate (leftS(ghostWidth,jsize-2*ghostWidth,nbVar))
   allocate (leftR(ghostWidth,jsize-2*ghostWidth,nbVar)) 
 
   allocate (rightS(ghostWidth,jsize-2*ghostWidth,nbVar))
   allocate (rightR(ghostWidth,jsize-2*ghostWidth,nbVar))

   end subroutine initS

   subroutine updateS(data)

   implicit none 

   integer :: i, j, iVar
   real(fp_kind), dimension(isize, jsize, nbVar), intent(in) :: data

   do iVar = 1, nbVar
       do i = 1,isize-2*ghostWidth
          do j = 1,ghostWidth
              topS(i,j,iVar) = data(i+ghostWidth, j+jsize-2*ghostWidth,iVar)
              bottomS(i,j,iVar) = data(i+ghostWidth, j+ghostWidth, iVar)
          end do
       end do

       do i = 1,ghostWidth
          do j = 1,jsize-2*ghostWidth
              leftS(i,j,iVar) = data(i+ghostWidth,j+ghostWidth , iVar)
              rightS(i,j,iVar) = data(i+isize-2*ghostWidth, j+ghostWidth, iVar)
          end do
         end do
   end do

  size_top = (isize-2*ghostwidth)*ghostWidth*nbVar
  size_right = (jsize-2*ghostwidth)*ghostWidth*nbVar

  end subroutine updateS

  subroutine comm

  integer :: nb_loop, index_loop,index_x,index_y

!write(*,*) size_x_max , " " , size_y_max

  nb_loop = size_y_max-1

!top send
  do index_loop = 0,nb_loop
   index_y = index_loop
   index_x = coord_x
   if (coord_y > 0 .AND. coord_y == index_y) then
    !write(*,*) index_x+size_x_max*(index_y),'receiving bottom from', index_x+size_x_max*(index_y-1)
    call MPI_RECV(bottomR, size_top, MPI_REAL, index_x+size_x_max*(index_y-1),1, MPI_COMM_WORLD, mpistat, ierr)
    !write(*,*) index_x+size_x_max*(index_y),'receiving bottom from', index_x+size_x_max*(index_y-1)
   end if 
   if (coord_y < nb_loop .AND. coord_y == index_y ) then
    call MPI_SEND(topS, size_top, MPI_REAL,index_x+size_x_max*(index_y+1) , 1, MPI_COMM_WORLD, mpistat, ierr)
   ! write(*,*) index_x+size_x_max*(index_y),'sending top to', index_x+size_x_max*(index_y+1) 
   end if
  end do

!write(*,*) 'envoi haut reussi'

!bottom send
  do index_loop = 0,nb_loop
  index_y = nb_loop - index_loop
  index_x = coord_x
   if (coord_y < nb_loop .AND. coord_y == index_y) then
    !write(*,*) coord_y,index_loop, index_x+size_x_max*(index_y),'receiving top from', index_x+size_x_max*(index_y+1) 
    call MPI_RECV(topR, size_top, MPI_REAL, index_x+size_x_max*(index_y+1), 2, MPI_COMM_WORLD, mpistat, ierr)
    !write(*,*) coord_y,index_loop, index_x+size_x_max*(index_y),'receiving top from', index_x+size_x_max*(index_y+1) 
   end if 
   if (coord_y > 0 .AND. coord_y == index_y ) then
    !write(*,*) index_x*size_y_max+coord_y ,'sending bottom to', (index_x-1)*size_y_max+coord_y
    call MPI_SEND(bottomS, size_top, MPI_REAL, index_x+size_x_max*(index_y-1), 2, MPI_COMM_WORLD, mpistat, ierr) 
    !write(*,*) index_x+size_x_max*index_y ,'sending bottom to', (index_y-1)*size_x_max+index_x
   end if
  end do

!write(*,*) 'envoi bas reussi'


  nb_loop = size_x_max-1 

!top send
  do index_loop = 0,nb_loop
   index_x = index_loop
   index_y = coord_y
   if (coord_x > 0 .AND. coord_x == index_x) then
    !write(*,*) index_x+size_x_max*(index_y),'receiving bottom from', index_x+size_x_max*(index_y-1)
    call MPI_RECV(leftR, size_top, MPI_REAL, (index_x-1)+size_x_max*index_y , 1, MPI_COMM_WORLD, mpistat, ierr)
    !write(*,*) index_x+size_x_max*(index_y),'receiving left from', (index_x-1)+size_x_max*index_y
   end if 
   if (coord_x < nb_loop .AND. coord_x == index_x ) then
    call MPI_SEND(rightS, size_top, MPI_REAL,(index_x+1)+size_x_max*index_y , 1, MPI_COMM_WORLD, mpistat, ierr)
    !write(*,*) index_x+size_x_max*(index_y),'sending right to', (index_x+1)+size_x_max*index_y 
   end if
  end do

!write(*,*) 'envoi droit reussi'

!bottom send
  do index_loop = 0,nb_loop
  index_x = nb_loop - index_loop
  index_y = coord_y
   if (coord_x < nb_loop .AND. coord_x == index_x) then
    !write(*,*) coord_y,index_loop, index_x+size_x_max*(index_y),'receiving top from', index_x+size_x_max*(index_y+1) 
    call MPI_RECV(rightR, size_top, MPI_REAL, (index_x+1)+size_x_max*index_y, 2, MPI_COMM_WORLD, mpistat, ierr)
    !write(*,*) index_x+size_x_max*(index_y),'receiving right from', (index_x+1)+size_x_max*index_y 
   end if 
   if (coord_x > 0 .AND. coord_x == index_x ) then
    !write(*,*) index_x*size_y_max+coord_y ,'sending bottom to', (index_x-1)*size_y_max+coord_y
    call MPI_SEND(leftS, size_top, MPI_REAL, (index_x-1)+size_x_max*index_y, 2, MPI_COMM_WORLD, mpistat, ierr) 
    !write(*,*) index_x+size_x_max*index_y ,'sending left to', (index_x-1)+size_x_max*index_y
   end if
  end do

!write(*,*) 'envoi gauche reussi'


   end subroutine comm

end module Communication
