module Partitioner

  use mpi

contains

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
      !write(*,*) 'The numbers of processors is equal to 2 to the power ', nbpowers

    end subroutine nbpow2

    subroutine partition(nbcores,nx,ny,mx,my,sizes_x,sizes_y)

      implicit none
      integer,intent(inout) :: nbcores, nx, ny, mx, my
      integer,intent(out),dimension(:),allocatable :: sizes_x, sizes_y
      integer :: nbpowers, modu, powx, powy, size_x, size_y, waste_x, waste_y, i

      call nbpow2(nbcores,nbpowers)
      modu = modulo(nbpowers,2)
      !write(*,*) 'modu is equal to',modu

      if ( modu == 1 ) then 
        powx = (nbpowers+1)/2
        powy = (nbpowers-1)/2
      else 
        powx = (nbpowers/2)
        powy = (nbpowers/2)
      end if

      mx = 2**powx
      my = 2**powy

      !write(*,*) 'mx is equal to', mx
      !write(*,*) 'my is equal to', my

      allocate (sizes_x(mx))
      allocate (sizes_y(my))

      waste_x = mod(nx,mx)
      size_x = (nx-waste_x)/mx
      waste_y = mod(ny,my)
      size_y = (ny-waste_y)/my

      do 10 i=1,mx
        if (waste_x > 0) then
          sizes_x(i) = size_x + 1
          waste_x = waste_x-1
        else 
          sizes_x(i) = size_x
        end if
        write(*,*) 'The size for x is', size_x
      10 continue

      do 20 i=1,my
        if (waste_y > 0) then
          sizes_y(i) = size_y + 1
          waste_y = waste_y-1
        else
          sizes_y(i) = size_y
        end if
        write(*,*) 'The size for y is', size_y
      20 continue

    end subroutine partition


end module Partitioner


