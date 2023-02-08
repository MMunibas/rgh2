(1.) 2A' potential energy surface for Ne---H2+ system.

Requirements:

RKHS toolkit - Download from https://github.com/MeuwlyGroup/RKHS

Compile a particular PES

A test program file (pes_test.f90) is given and can be compiled as

gfortran RKHS.f90 Ne--H2+_PV5Z_CCSDT.f90 pes_test.f90

Running the executable

Before running the executable make sure that the .coeff, .csv and/or .kernel files for that PES present in the current directory (or change the file path in the fortran program).
