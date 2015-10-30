# ######## GNU ########
 F90 = gfortran
 MPIF90 = mpif90
 CFLAGS = -C -g -Wall
 FFLAGS = -O3 -Wall
 #FFLAGS = -g -Wall
 # -ffpe-trap=zero,invalid,underflow -fbacktrace
 LDFLAGS = -O3

# ######## INTEL ########
# F90 = ifort
# FFLAGS = -warn all -O3 -xSSE4.2 -ipo
# #FFLAGS = -xSSE4.2 -ipo 
# #FFLAGS = -O3 -no-vec
# LDFLAGS = -O3

# ######## PGI ########
#F90 = pgfortran
#FFLAGS = -O4 -fast -Minform=warn -Minfo=all -Mpreprocess 
##-Mipa=inline -Mvect=simd:256
##FFLAGS = -g -Minform=warn -Minfo -Mpreprocess -Ktrap=fp
#LDFLAGS = -O4
##LDFLAGS = -g

SRCDIR = .
OBJ = \
	m_HydroPrecision.o \
	m_HydroConstants.o \
	m_HydroParameters.o \
	m_HydroUtils.o \
	m_Monitoring.o \
	m_HydroRun.o \
	m_Partitioner.o \
	main.o

euler2d: $(OBJ)
#$(F90) $(LDFLAGS) $(OBJ) -o $@
	 $(MPIF90) $(CFLAGS) $(OBJ) -o $@

clean:
	rm -f *.o *.mod euler2d

cleanall: clean
	rm -f *.vti

%.o:    $(SRCDIR)/%.f90
	$(MPIF90) $(FFLAGS) -c $<

