#!/bin/bash

#BSUB -q long
#BSUB -W 200:00
#BSUB -R rusage[mem=12000]
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -e output.err
#BSUB -oo output.log

module load R/3.6.1
module load gcc/8.1.0


Rscript my_RAD_parse_test.R all_urpa.mafs
