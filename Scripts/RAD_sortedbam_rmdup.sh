#!/bin/bash
#Remove PCR duplicates from sorted bam files
#BSUB -q short
#BSUB -W 4:00
#BSUB -R rusage[mem=1600]
#BSUB -n 6
#BSUB -R span[hosts=1]
#BSUB -e rapture_rmdups.err
#BSUB -oo rapture_rmdups.log


module load samtools/1.9
mkdir ../bam_sort_fltr_wA02
for file in *_sort.bam

do

sample=`echo $file | cut -f1,2,3,4 -d "_" `

echo $sample

samtools rmdup "$sample"_sort.bam ../bam_sort_fltr_wA02/"$sample"_sortfltr.bam

done
