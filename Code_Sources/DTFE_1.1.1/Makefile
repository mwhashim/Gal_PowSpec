# Makefile for compiling the DTFE code on Linux systems


# Path to the GSL, Boost C++ and CGAL libraries - must be set by user (only if they aren't installed in the default system path) -- (NOTE: You must add only the directory where the libraries are installed, the program will add the '/lib' and '/include' parts automatically); C++ compiler - preferably a version that supports OpenMP
GSL_PATH = /opt/apps/libs/gsl/2.1/gcc-4.4.7 
BOOST_PATH = /opt/gridware/pkg/libs/boost/1.53.0/gcc-4.4.7+openmpi-1.8.1+python-2.7.5  
CGAL_PATH = /opt/apps/cgal/4.8   
CC = g++


# paths to where to put the object files and the executables files. If you build the DTFE library than you also need to specify the directory where to put the library and the directory where to copy the header files needed by the library (choose an empty directory for the header files).
OBJ_DIR = ./o
BIN_DIR = ./
LIB_DIR = ./
INC_DIR = ./DTFE_include


# Choose the compiler directives
OPTIONS = 
#------------------------ set the number of spatial dimensions (2 or 3 dimensions)
OPTIONS += -DNO_DIM=3 
#------------------------ set type of variables - float (comment the next line) or double (uncomment the next line)
# OPTIONS += -DDOUBLE 
#------------------------ set which quantities can be computed (can save memory by leaving some out)
# Comment this line if you don't need to compute velocity and velocity related components 
OPTIONS += -DVELOCITY 
# Comment this line if you don't need to interpolate additional fields stored in the scalar variable
# OPTIONS += -DSCALAR 
# number of components of the scalar variable
# OPTIONS += -DNO_SCALARS=3 
#------------------------ additional compiler options
# enable to check if the padding gives a complete Delaunay Tesselation of the region of interest
# OPTIONS += -DTEST_PADDING 
# enable this option if to use OpenMP (share the workload between CPU cores sharing the same RAM)
OPTIONS += -DOPEN_MP 
# enable this if the input data is not in Mpc units (value=what is 1 MPc in the units of the data - in this example the data is in kpc)
# OPTIONS += -DMPC_UNIT=1000. 
#------------------------ options usefull when using DTFE as a library
# uncomment the line to get access to a function that returns the Delaunay triangulation of the point set
# OPTIONS += -DTRIANGULATION 
#-------- special customizations for Wojtek - a dedicated user of DTFE
#OPTIONS += -DWOJTEK 


#------------------------ compiler directive that affect only the help messages when using the '-h / --help' option (it does not affect the program in any other way)- if the option is uncommented, than it will show that set of options
OPTIONS += -DFIELD_OPTIONS 
OPTIONS += -DREGION_OPTIONS 
# OPTIONS += -DPARTITION_OPTIONS 
# OPTIONS += -DPADDING_OPTIONS 
OPTIONS += -DAVERAGING_OPTIONS 
# OPTIONS += -DREDSHIFT_CONE_OPTIONS 
OPTIONS += -DADDITIONAL_OPTIONS 







###############  DO NOT MODIFY BELOW THIS LINE  ###########################
# do not modify below this line
SRC = ./src
INCLUDES = 
LIBRARIES = 

ifneq ($(strip $(GSL_PATH)),)
	INCLUDES += -I/$(strip $(GSL_PATH))/include 
	LIBRARIES += -L/$(strip $(GSL_PATH))/lib 
endif
ifneq ($(strip $(BOOST_PATH)),)
	INCLUDES += -I/$(strip $(BOOST_PATH))/include 
	LIBRARIES += -L/$(strip $(BOOST_PATH))/lib 
endif
ifneq ($(strip $(CGAL_PATH)),)
	INCLUDES += -I/$(strip $(CGAL_PATH))/include 
	LIBRARIES += -L/$(strip $(CGAL_PATH))/lib 
endif

COMPILE_FLAGS = -frounding-math -O3 -fopenmp -DNDEBUG $(OPTIONS)
DTFE_INC = $(INCLUDES)
# the following libraries should work in most cases
DTFE_LIB = -rdynamic $(LIBRARIES) -lCGAL -lboost_thread -lboost_filesystem -lboost_program_options -lgsl -lgslcblas -lm -lboost_system
# the following libraries work on "Fedora 15"
# DTFE_LIB = -rdynamic $(LIBRARIES) -lCGAL -lboost_thread-mt -lboost_filesystem -lboost_program_options -lgsl -lgslcblas -lgmp -lboost_system



IO_SOURCES = $(addprefix io/, input_output.h gadget_reader.cc text_io.cc binary_io.cc my_io.cc)
MAIN_SOURCES = main.cpp DTFE.h message.h user_options.h input_output.cc $(IO_SOURCES)
DTFE_SOURCES = DTFE.cpp define.h particle_data.h user_options.h box.h quantities.h user_options.cc quantities.cc subpartition.h random.cc CIC_interpolation.cc TSC_interpolation.cc SPH_interpolation.cc kdtree/kdtree2.hpp Pvector.h message.h
TRIANG_SOURCES = $(addprefix CGAL_triangulation/, triangulation.cpp triangulation_miscellaneous.cc unaveraged_interpolation.cc averaged_interpolation_1.cc averaged_interpolation_2.cc padding_test.cc CGAL_include_2D.h CGAL_include_3D.h vertexData.h particle_data_traits.h) define.h particle_data.h user_options.h box.h quantities.h Pvector.h message.h math_functions.h

ALL_FILES = $(DTFE_SOURCES) $(TRIANG_SOURCES) $(MAIN_SOURCES) kdtree/kdtree2.hpp kdtree/kdtree2.cpp
LIB_FILES = $(DTFE_SOURCES) $(TRIANG_SOURCES)

HEADERS_1 = DTFE.h define.h user_options.h particle_data.h quantities.h Pvector.h math_functions.h 
HEADERS_2 = $(addprefix CGAL_triangulation/, CGAL_include_2D.h CGAL_include_3D.h vertexData.h particle_data_traits.h)



DTFE: set_directories $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation.o $(OBJ_DIR)/main.o $(OBJ_DIR)/kdtree2.o Makefile
	$(CC) $(COMPILE_FLAGS) $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation.o $(OBJ_DIR)/main.o $(OBJ_DIR)/kdtree2.o $(DTFE_LIB) -o $(BIN_DIR)/DTFE


$(OBJ_DIR)/main.o: $(addprefix $(SRC)/, $(MAIN_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/main.o -c $(SRC)/main.cpp

$(OBJ_DIR)/DTFE.o: $(addprefix $(SRC)/, $(DTFE_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/DTFE.o -c $(SRC)/DTFE.cpp

$(OBJ_DIR)/kdtree2.o: $(SRC)/kdtree/kdtree2.hpp $(SRC)/kdtree/kdtree2.cpp Makefile
	$(CC) -O3 -ffast-math -fomit-frame-pointer $(DTFE_INC) -o $(OBJ_DIR)/kdtree2.o -c $(SRC)/kdtree/kdtree2.cpp

$(OBJ_DIR)/triangulation.o: $(addprefix $(SRC)/, $(TRIANG_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/triangulation.o -c $(SRC)/CGAL_triangulation/triangulation.cpp


library: set_directories set_directories_2 $(addprefix $(SRC)/, $(LIB_FILES) ) copy_headers Makefile
	$(CC) $(COMPILE_FLAGS) -fPIC $(DTFE_INC) -o $(OBJ_DIR)/DTFE_l.o -c $(SRC)/DTFE.cpp
	$(CC) -O3 -ffast-math -fomit-frame-pointer -fPIC $(DTFE_INC) -o $(OBJ_DIR)/kdtree2_l.o -c $(SRC)/kdtree/kdtree2.cpp
	$(CC) $(COMPILE_FLAGS) -fPIC $(DTFE_INC) -o $(OBJ_DIR)/triangulation_l.o -c $(SRC)/CGAL_triangulation/triangulation.cpp
	$(CC) $(COMPILE_FLAGS) -shared $(OBJ_DIR)/DTFE_l.o $(OBJ_DIR)/triangulation_l.o $(OBJ_DIR)/kdtree2_l.o $(DTFE_LIB) -o $(LIB_DIR)/libDTFE.so


clean:
	rm -f $(BIN_DIR)/DTFE $(OBJ_DIR)/*.o

copy_headers:
	cp $(addprefix $(SRC)/, $(HEADERS_1)) $(INC_DIR)
	cp $(addprefix $(SRC)/, $(HEADERS_2)) $(INC_DIR)/CGAL_triangulation

set_directories:
	@ if !( test -d $(OBJ_DIR) ); \
	then mkdir $(OBJ_DIR); \
	fi
	@ if !( test -d $(BIN_DIR) ); \
	then mkdir $(BIN_DIR); \
	fi

set_directories_2:
	@ if !( test -d $(LIB_DIR) ); \
	then mkdir $(LIB_DIR); \
	fi
	@ if !( test -d $(INC_DIR) ); \
	then mkdir $(INC_DIR); \
	fi
	@ if !( test -d $(INC_DIR)/CGAL_triangulation ); \
	then mkdir $(INC_DIR)/CGAL_triangulation; \
	fi
