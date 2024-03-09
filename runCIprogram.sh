#!/bin/bash

gfortran -c functions/General.f95
gfortran -c functions/Configuration.f95
gfortran -c functions/LoadInput.f95
gfortran -c functions/Hamiltonian.f95
gfortran -c functions/CI.f95

gfortran -o prog CIprogram.f95 CI.o hamiltonian.o loadInput.o configuration.o general.o -llapack
./prog