#!/bin/bash


cat > ilqins_4x2x12.key <<EOF
parameters amoebabio18.prm
verbose

a-axis 120.0
b-axis 120.0
c-axis 120.0

openmp-threads 4           # No. of parallel CPU threads to use

cutoff 9.0                 # Nonbonded interactions direct cutoff
ewald                      # Switch on PME
ewald-cutoff 7.0           # PME real space cutoff (overrides the 9A above)
vdw-correction             # Switch on analytical long-range vdW correction

polar-eps 0.01             # Dipole convergence criterion in RMS Debye/atm

neighbor-list              # Use pairwise neighbor list for calculating
                           # nonbonded interactions (improves speed)
EOF

##this command takes the xyz file we already built
##and gives an info dump that can be a starting point to
##build a special amber-amoeba parmtop file.
echo "running tinker analyze"
/home/josh/tinker/bin/analyze -k ilqins_4x2x12.key ilqins_4x2x12_ipqVac.xyz PC > ilqins_4x2x12.analout

echo "build analout: "
ls -l *analout


cat > ilqins_4x2x12.key <<EOF
parameters amoebabio18
verbose
a-axis 120.0
b-axis 120.0
c-axis 120.0
cutoff 9.0                 # Nonbonded interactions direct cutoff
ewald                      # Switch on PME
ewald-cutoff 7.0           # PME real space cutoff (overrides the 9A above)
vdw-correction             # Switch on analytical long-range vdW correction

polar-eps 0.01             # Dipole convergence criterion in RMS Debye/atm

neighbor-list              # Use pairwise neighbor list for calculating
                           # nonbonded interactions (improves speed)
EOF


##tinker_to_amber cannot overwrite files, so make sure to clear old versions:
rm -f ilqins_4x2x12_tinker2amber.top
rm -f ilqins_4x2x12_tinker2amber.crd

##build the amber inputs, a special amoeba topfile and an associated
##crd.
$AMBERHOME/bin/tinker_to_amber  -key ilqins_4x2x12.key <<EOF
ilqins_4x2x12.analout
ilqins_4x2x12.xyz
ilqins_4x2x12.pdb
ilqins_4x2x12
ilqins_4x2x12_tinker2amber.top
ilqins_4x2x12_tinker2amber.crd
EOF
			      

