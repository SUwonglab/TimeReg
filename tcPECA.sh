#/bin/bash

# Time Course Regulatory Analysis(TCRA) updated Aug 14th 2019

cat TCRA_config TCRA.m > run_TCRA.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "run_TCRA; exit"
