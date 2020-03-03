#/bin/bash

# Time Course Regulatory Analysis(TCRA) updated Aug 14th 2019
module load matlab
cat TCRA_config TCRA.m > run_TCRA.m
matlab -nodisplay -nosplash -nodesktop -r "run_TCRA; exit"
