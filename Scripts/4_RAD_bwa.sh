#!/bin/bash
#run script in directory where files are, or change path accordingly below; would need to also change cut command accordingly
#BSUB -q long
#BSUB -W 40:00
#BSUB -R rusage[mem=8000]
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -e rad_bwa.err
#BSUB -oo rad_bwa.log

module load bwa/0.7.5a
module load samtools/1.9

REF=project/uma_lisa_komoroske/path/to/ref #change this to be path to your reference genome file

bwa index $REF

mkdir ./samfiles

for file in *_RA.fastq

do
echo $file

sample=`echo $file |cut -f1,2,3,4,5 -d "_"` #change the cut fields depending on how many fields are in your sample name, for example my files have 5 identifiers separated by '_' before the RA.fastq ending

echo $sample #confirm that your sample name contains all sample name fields that you want to keep in the next step, not file ending

bwa mem -t 20 $REF "$sample"_RA.fastq  "$sample"_RB.fastq | samtools sort -O bam -o ../4_mapped_bam_files/"$sample"_sort.bam \
2> ../4_mapped_bam_files/"$sample"_sort.stderr

done



mkdir ../4_mapped_bam_files/sort_flagstat
cd ../4_mapped_bam_files/
for file in *.bam

do

sample=`echo $file | cut -f1 -d "." `

echo $sample
samtools flagstat "$sample".bam \
> ./sort_flagstat/"$sample"_flagstat.txt \
2> ./sort_flagstat/"$sample"_flgstbam.stderr

done


#this is just combining all the bam stats files for each sample into one text file that is tab delimited so one sample per row so can then easily manipulate to summarize mapping stats etc

cd ./sort_flagstat/
for filename in *_flagstat.txt; do
    cat "$filename"|tr '\n' '\t'
    echo "$filename"
done > ./All_RAD_sort_combined.flagstat_tabbed.txt


#few lines to generate a cleaned up summary file of mapping stats for all samples to compare pre and post filtering

awk '{print NF}' All_RAD_sort_combined.flagstat_tabbed.txt #88 columns
awk '{print $88 "\t" $1 "\t" $23 "\t" $27 "\t"$30 "\t"$44 "\t" $49 "\t" $52 "\t"$60 "\t" $64 "\t" $67 "\t" $77 }' All_RAD_sort_combined.flagstat_tabbed.txt > All_RAD_sort_combined_flagstat_reformat.txt
sed -i 's/(//g' All_RAD_sort_combined_flagstat_reformat.txt #-i changes it in the original file; could test sending to another outfile first
sed -i 's/%//g' All_RAD_sort_combined_flagstat_reformat.txt
cat ../../scripts_and_keyfiles/short_flagstat_headers All_RAD_sort_combined_flagstat_reformat.txt >All_head_flagstat_reformat.txt
