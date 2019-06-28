FC = mpifort  # or mpiifort

HDF5_DIR     = /home/ohull/mylibraries/gnu
HDF5_LDIR    =  $(HDF5_DIR)/lib
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
								$(HDF5_LDIR)/libhdf5_hl.a \
							  $(HDF5_LDIR)/libhdf5_fortran.a \
							  $(HDF5_LDIR)/libhdf5.a -lz -ldl

HDF5INCLUDE  = $(HDF5_DIR)/include

FCFLAGS = -I $(HDF5INCLUDE)
LDFLAGS = $(HDF5LIB)

all: my_main.x

%.x: %.f90
		$(FC) -g $^ -o $@ $(FCFLAGS) $(LDFLAGS)

clean:
		rm -rf *.x *.o

.PHONY: clean
