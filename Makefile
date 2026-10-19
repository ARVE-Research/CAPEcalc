# makefile

# 1. specify the name of your Fortran compiler

FC = gfortran

# 2. specify compiler options, including openMP parallelization if desired

# optimized flags
FCFLAGS = -O3 -ffree-form -ffree-line-length-none -ftree-vectorize -fopenmp -Wall

# debugging flags
# FCFLAGS = -g -Og -ffree-line-length-none -fcheck=all -fno-check-array-temporaries -ffpe-trap=invalid,zero,overflow,underflow -g -fbacktrace -Wall -pedantic

# 3. provide the path to your netCDF-Fortran libraries

netcdf = /work/kaplan_lab/easybuild/software/netCDF-Fortran/4.6.1-gompi-2024a

# should not need to modify anything below this line

# ---------------------------------------------

NC_LIB = $(netcdf)/lib
NC_INC = $(netcdf)/include

CPPFLAGS = -I$(NC_INC)
LDFLAGS  = -L$(NC_LIB)
LIBS     = -lnetcdff

# ---------------------------------------------

CAPE_OBJS = parametersmod.o \
            capemod.o       \
            calendarmod.o   \
            CAPEcalc.o

# ---------------------------------------------

.SUFFIXES: .o .f90 .F90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.F90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	CAPEcalc

CAPEcalc: $(CAPE_OBJS)
	$(FC) $(FCFLAGS) -o CAPEcalc $(CAPE_OBJS) $(LDFLAGS) $(LIBS)

clean::	
	-rm *.o *.mod CAPEcalc
