#!/bin/bash

MOD=$(pwd)/modules
SRC=$(pwd)/src
OBJ=$(pwd)/Output_Files

USRMOD=$(pwd)/user_codes/modules
USRSRC=$(pwd)/user_codes/src

rm -r $OBJ/*.*

gfortran -J$OBJ -c $MOD/Types.f90 -o $OBJ/Types.o
gfortran -J$OBJ -c $MOD/ParamIO.f90 -o $OBJ/ParamIO.o
gfortran -J$OBJ -c $MOD/Linkedlist_Handling.f90 -o $OBJ/Linkedlist_Handling.o
gfortran -J$OBJ -c $MOD/Bandwidth.f90 -o $OBJ/Bandwidth.o
gfortran -J$OBJ -c $MOD/Boundaryconditions.f90 -o $OBJ/Boundaryconditions.o
gfortran -J$OBJ -c $MOD/Controlparameters.f90 -o $OBJ/Controlparameters.o
gfortran -J$OBJ -c $MOD/Dynamicstepparameters.f90 -o $OBJ/Dynamicstepparameters.o
gfortran -J$OBJ -c $MOD/Globals.f90 -o $OBJ/Globals.o
gfortran -J$OBJ -c $MOD/Mesh.f90 -o $OBJ/Mesh.o
gfortran -J$OBJ -c $MOD/Staticstepparameters.f90 -o $OBJ/Staticstepparameters.o
gfortran -J$OBJ -c $MOD/Printparameters.f90 -o $OBJ/Printparameters.o
gfortran -J$OBJ -c $MOD/Fileparameters.f90 -o $OBJ/Fileparameters.o
gfortran -J$OBJ -c $MOD/Stiffness.f90 -o $OBJ/Stiffness.o
gfortran -J$OBJ -c $MOD/User_Subroutine_Storage.f90 -o $OBJ/User_Subroutine_Storage.o

gfortran -J$OBJ -c $USRMOD/Element_Utilities.f90 -o $OBJ/Element_Utilities.o

cp $USRSRC/aba_param.inc $OBJ
cp $USRSRC/vaba_param.inc $OBJ

cd $OBJ

gfortran -c $USRSRC/*.f90 -I$OBJ
gfortran -c $USRSRC/*.for -I$OBJ

gfortran -c $SRC/*.f90 -I$OBJ

gfortran -o main *.o