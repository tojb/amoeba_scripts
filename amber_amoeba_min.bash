#!/bin/bash

cat > amoebaHeat.in <<EOF
test of amoeba
&cntrl
   iamoeba=1,
   ntx=1,
   irest=0,
   dt=0.0005, nstlim=2000, 
   ntpr=1, ntwx=20,
   vlimit=20,
   cut=8.0, maxcyc=100,
   ntt=3, tempi=0., temp0=300., gamma_ln=10.0,
/
&ewald
   order=6,
/
&amoeba
   dipole_scf_tol = 0.01, dipole_scf_iter_max=30,
/
EOF


step=amoebaHeat
sys=ilqins_4x2x12_tinker2amber

#mpirun -np 4 $AMBERHOME/bin/pmemd.amoeba.MPI \

$AMBERHOME/bin/pmemd.amoeba \
       -O -i $step.in -p $sys.top -c $sys.crd \
         -o  $step.out -x $step.x -r $step.rst &
