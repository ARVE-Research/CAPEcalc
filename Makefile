# makefile

FC = gfortran
# FCFLAGS = -ffree-form -ffree-line-length-none -ftree-vectorize -Wall

FCFLAGS = -g  -O0 -ffree-line-length-none -fcheck=all -fno-check-array-temporaries -ffpe-trap=invalid,zero,overflow,underflow -g -fbacktrace -Wall -pedantic

# use the command "nf-config --all" to find the location of your netCDF installation
# and enter the path next to " --prefix    ->" on the line below

netcdf = /usr/local

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
            merra2cape.o

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

all::	merra2cape

merra2cape: $(CAPE_OBJS)
	$(FC) $(FCFLAGS) -o merra2cape $(CAPE_OBJS) $(LDFLAGS) $(LIBS)

clean::	
	-rm *.o *.mod merra2cape
