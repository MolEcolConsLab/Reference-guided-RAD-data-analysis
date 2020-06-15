#!/bin/bash
#run script in directory where files are, or change path accordingly below; would need to also change cut command accordingly
#BSUB -q long
#BSUB -W 40:00
#BSUB -R rusage[mem=8000]
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -e rad_bwa.err
#BSUB -oo rad_bwa.log

module load bwa/0.7.
module load samtools/1.9
bwa index /project/uma_lisa_komoroske/GT_genome/Cmyd.v1.1.fa
mkdir ./samfiles
for file in ./*_RA.fastq

do
echo $file

sample=`echo $file |cut -f1,2,3,4 -d "_"` #change the cut fields depending on how many fields are in your sample name

echo $sample #confirm that your sample name contains all sample name fields that you want to keep in the next step, not file ending

bwa mem -t 20 /project/uma_lisa_komoroske/GT_genome/Cmyd.v1.1.fa "$sample"_RA.fastq  "$sample"_RB.fastq | samtools sort -O bam -o ../4_mapped_bam_files/"$sample"_sort.bam \
2> ../4_mapped_bam_files/"$sample"_sort.stderr

done

wait


mkdir ../4_mapped_bam_files/sort_flagstat
for file in *.bam

do

sample=`echo $file | cut -f1 -d "." `

echo $sample
samtools flagstat "$sample".bam \
> ../4_mapped_bam_files/sort_flagstat/"$sample"_flagstat.txt \
2> ../4_mapped_bam_files/sort_flagstat/"$sample"_flgstbam.stderr

done
