!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!< Defines all variables needed to perform parallel partitioning
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Partitioner

  use mpi

  ! partitioner parameters
  integer,dimension(:),allocatable :: sizes_x, sizes_y

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! for 2^n processors, determine n
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine nbpow2(nbcores,nbpowers)

  implicit none

      integer,intent(inout) :: nbcores
      integer,intent(out) :: nbpowers

      nbpowers = 1
      do  
        nbcores = nbcores/2
        if (nbcores <= 1) exit
        nbpowers = nbpowers + 1
      end do

    end subroutine nbpow2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! do the partitioning
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine partition(nbcores,nx,ny)

      implicit none

      ! dummy variables
      integer,intent(inout) :: nbcores, nx, ny
      ! local variables
      integer :: nbpowers, modu, powx, powy, mx, my, size_x, size_y, waste_x, waste_y, i

      ! determine how many cells in each direction
      call nbpow2(nbcores,nbpowers)
      modu = modulo(nbpowers,2)

      if ( modu == 1 ) then 
        powx = (nbpowers+1)/2
        powy = (nbpowers-1)/2
      else 
        powx = (nbpowers/2)
        powy = (nbpowers/2)
      end if

      mx = 2**powx
      my = 2**powy

      ! in case there are extra cells
      waste_x = mod(nx,mx)
      size_x = (nx-waste_x)/mx
      waste_y = mod(ny,my)
      size_y = (ny-waste_y)/my

      allocate (sizes_x(mx))
      allocate (sizes_y(my))

      do i=1,mx
        if (waste_x > 0) then
          sizes_x(i) = size_x + 1
          waste_x = waste_x-1
        else 
          sizes_x(i) = size_x
        end if
      end do

      do i=1,my
        if (waste_y > 0) then
          sizes_y(i) = size_y + 1
          waste_y = waste_y-1
        else
          sizes_y(i) = size_y
        end if
      end do

    end subroutine partition


end module Partitioner


