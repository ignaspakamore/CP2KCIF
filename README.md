# CP2KCIF
CP2K input/restart file converter to cif or xyz

out2cif automaticly detects whether coordinates are scaled or not and converts cartesian to fraction or vice versa for cif and xyz files.
Works only for units in angstrom. 

use: python3 out2cif.py -cif/-xyz cp2k.inp

