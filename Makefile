#fopts= -fdefault-real-8 -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid
export incl=-I yasso/fortran
fopts= -fdefault-real-8 -fPIC -O3 $(incl)
#export F90C=gfortran $(fopts)

export SRCDIR=model/
export OBJDIR=$(SRCDIR)
export PREPROCESS=-x f95-cpp-input
export FFLAGS=-I $(OBJDIR)
export F90C=gfortran $(fopts) $(FFLAGS)

export srcfiles=$(shell find model -name *.f90)

# Compile basgra into a shared library callable from R or Fortran.
basgra : dep yasso
	$(MAKE) -f Makefile.basgra BASGRA.DLL

# Compile basgra into a python module using f2py
basgrafort.so : dep yasso
	$(MAKE) -f Makefile.basgra $@

.PHONY: yasso
# Build Yasso into a static library
yasso :
	$(MAKE) -C yasso libyasso.a

.PHONY: dep clean
# Find the dependencies within BASGRA
dep : 
	find model -name *.f90 > filelist
	perl config.pl filelist BASGRA null no_mpi

clean :
	rm -f model/*.o model/*.mod

