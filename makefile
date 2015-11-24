PETSC_DIR := /usr/local/lib/petsc/3.6.1
	
CFLAGS	        =
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =
EXAMPLESC       = 


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

all: ex12

ex12: ex12.o   chkopts
	-${CLINKER} -o ex12 ex12.o  ${PETSC_DM_LIB}
	${RM} -f ex12.o



