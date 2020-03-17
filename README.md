# amoeba_scripts
Some attempts to run amoeba using tinker / amber.

Workflow is:

## pdbxyz_amber_jtb.py

Use this script to convert an amber-format PDB file to an amoeba xyz.  Works with many chains (does not read chain ID from the PDB), which not all pdbxyz executables seem to do.

## build_amber.bash

It turns out that current (2020) downloads of tinker/analyze and amber/tinker_to_amber are OK as far as I can see, so this bash script is just a record of using these tools to dump out forcefield info (using analyze) and repack that into a non-standard amber parmtop format for pmemd.amoeba to use (tinker_to_amber).

## amber_amoeba_min.bash

This is an MD run with a small timestep and langevin thermostat, basically a test / proof of concept establishing that pmemd.amoeba is working for my system.  I was prompted by a warning to increase the spline order, I'm not sure if this is general to amoeba runs or just a quirk of my system.

I've verified also that outputs from this toolchain work with tinker's minimize program.





