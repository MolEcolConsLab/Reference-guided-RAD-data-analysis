#!/bin/bash
#run fastqc

#BSUB -q short
#BSUB -W 4:00 
#BSUB -R rusage[mem=1000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e fastqcRA.err
#BSUB -oo fastqcRA.log

module load fastqc/0.11.5

for file in *RA.fastq.gz

do

echo $file

fastqc $file

done

