PETSC_DIR := /usr/local/lib/petsc/3.6.1
	
CFLAGS	        = 
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =
EXAMPLESC       = 


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

all: seismic

seismic: seismic.o vti_io.o  chkopts
	-${CLINKER} -o seismic seismic.o vti_io.o ${PETSC_DM_LIB}
	${RM} -f *.o



