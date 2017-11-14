#!/bin/bash --login

cd ~/Downloads/plumed-test/
make clean
./configure 
make -j 8
source sourceme.sh
cd ~/Downloads/lammps-2Apr13-replica/
plumed patch -r 
plumed patch -p
cd src/
rm Obj_openmpi/fix_plumed.o
make openmpi -j 8

