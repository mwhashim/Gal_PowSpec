#
# The compiler
# ------------
#
FC = ifort
#
# The key words of the compiler (optimization, etc)
# -----------------------------
#
# Activate the precompilation parameter OMP for calculations
# using openMP (-DOMP in most fortran90 compilers)
#
# Activate the precompilation parameter -DADD1US for
# fixing underscore problems while linking with fftw2.
#
# Activate big-endian/little-endian conversion compiling option when
# relevant to avoid problems while reading input simulations files.
#
FCFLAGS=  -Vaxlib -CB -mcmodel=large -O -openmp -DOMP
#
# The path of FFTW2 (NOT 3) 
# -----------------
# 
# Use -lfftw -lrfftw for standard calculation
#
# Use  -lrfftw_threads -lrfftw -lfftw_threads -lfftw -lpthread for
# calculations with openMP.
#
LIBS = -lrfftw_threads -lrfftw -lfftw_threads -lfftw -lpthread

# In principle, no need to change anything from here.
INC =

OBJS =  ../src/twopidef.o \
	../src/numrec_tools.o \
        ../src/random_number_tools.o \
        ../src/fourier_tools.o \
        ../src/gadget_ramses_tools.o \
        ../src/fourier_taylor_tools.o \
        ../src/fourier_tools3D.o \
        ../src/powmes_common.o \
        ../src/powmes.o

all:   powmes

powmes:  $(OBJS)
	 $(FC) $(FCFLAGS)  $(OBJS) -o $@ $(INC) $(LIBS)

clean: 
	/bin/rm  ../src/*.o *.mod 

.SUFFIXES:
.SUFFIXES: .f90 .F90 .o
.f90.o:
	$(FC) $(FCFLAGS)  $(INC) -c $< -o $@
.F90.o:
	$(FC) $(FCFLAGS)  $(INC) -c $< -o $@
