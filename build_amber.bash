#!/bin/bash

if [ "$#" -ge "1" ]
then
    sysName=$1
else
    echo "$0 : require system name (eg ilqins if reading ilqins.xyz)"
    exit 1
fi

tink="$HOME/tinker"
parmdir=$tink/params


cat > $sysName.key <<EOF
parameters $parmdir/amoebabio18.prm
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
$tink/bin/analyze -k $sysName.key $sysName.xyz PC > $sysName.analout

echo "build analout: "
ls -l *analout


cat > $sysName.key <<EOF
parameters $parmdir/amoebabio18
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

##build the amber inputs, a special amoeba topfile and an associated
##crd.

##tinker_to_amber cannot overwrite files, so make sure to clear old versions:
rm -f ${sysName}_tinker2amber.top
rm -f ${sysName}_tinker2amber.crd
$AMBERHOME/bin/tinker_to_amber  -key $sysName.key <<EOF
$sysName.analout
$sysName.xyz
$sysName.pdb
$sysName
${sysName}_tinker2amber.top
${sysName}_tinker2amber.crd
EOF
			      
echo "built: ${sysName}_tinker2amber.top  ${sysName}_tinker2amber.crd"
