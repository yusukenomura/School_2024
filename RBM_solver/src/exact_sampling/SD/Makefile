# Copyright (c) 2024 Yusuke Nomura

#======= please edit this part ========
FC = gfortran
FFLAGS =
LIBBLAS =
#======================================

.SUFFIXES: .f90
.f90.o:
	${FC} ${FFLAGS} -c $<

TARGET = RBM_solver_SD.x
OBJECTS = main.o sub.o mod_RBM.o
MOD_FILES = mod_rbm.mod


${TARGET} : ${OBJECTS} 
	${FC} ${FFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 

sub.o: sub.f90 mod_RBM.o 
main.o: main.f90 mod_RBM.o


clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} 
