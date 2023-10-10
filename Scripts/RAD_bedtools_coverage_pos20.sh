#!/bin/bash
# run bedtools coverage for file in directory *.bam
#SBATCH --time=7:00:00
#SBATCH -c 1
#SBATCH --mem=5G

module load bedtools/2


bamdir=/nese/meclab/Blair/Chemyd_Brazil/alignments/rmdupe

for file in $bamdir/*.bam

do

sample="$(basename -- $file .bam)" \

echo $sample
bedtools coverage -hist -a RAD_bedtools_pos20_min134_brazil.bed -b $file >./bedtools_pos20_cov/"$sample".pos20.cov.txt

done
