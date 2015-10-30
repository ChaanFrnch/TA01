module Partitioner

	use mpi

	contains

		subroutine nbpow2(nbcores,nbpowers)

			implicit none
			int,intent(in) :: nbcores
			int,intent(out) :: nbpowers

			nbpowers = 1
			do 
				nbcores = nbcores/2
				if (nbcores =< 1) exit
				nbpowers = nbpowers + 1
			end do

		end subroutine nbpow2

		subroutine partition(nbcores,nx,ny,mx,my,sizes_x,sizes_y)

			implicit none
			integer,intent(inout) 												   :: nbcores,nx,ny,mx,my
			integer,intent(out),dimension(:),allocatable :: sizes_x,sizes_y
			
			integer :: nbpowers,mod,powx,powy,_size_x,_size_y,waste_x,waste_y
			call nbpow2(nbcores,nbpowers)
			mod = modulo(nbpowers,2)

			if ( mod = 1 ) then 
				powx = (nbpowers+1)/2
				powy = (nbpowers-1)/2			
			else 
				powx = (nbpowers/2)
				powy = (nbpowers/2)			
			end if
	
			mx = 2**powx
			my = 2**powy

			call allocate(sizes_x(mx))
			call allocate(sizes_y(my))

			waste_x = mod(nx,mx)
			_size_x = (nx-waste_x)/mx
			waste_y = mod(ny,my)
			_size_y = (ny-waste_y)/my

			do 10 i=1,mx
				if (waste_x > 0) then
					sizes_x(i) = _size_x + 1
					waste = waste-1
				else 
					sizes_x(i) = _size_x
				end if
			10 continue

			do 10 i=1,my
				if (waste_x > 0) then
					sizes_y(i) = _size_x + 1
					waste = waste-1
				else
					sizes_x(i) = _size_x
				end if
			10 continue
				

		end subroutine partition


end module Partitioner

 
