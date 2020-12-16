#!/bin/bash
#Remove PCR duplicates from sorted bam files, then generate flagstats
#run from within directory containing sorted bam files
#BSUB -q short
#BSUB -W 4:00
#BSUB -R rusage[mem=1600]
#BSUB -n 6
#BSUB -R span[hosts=1]
#BSUB -e RAD_rmdups_mod.err
#BSUB -oo RAD_rmdups_mod.log



for file in *_sort.bam

do

sample=`echo $file | cut -f1,2,3,4,5 -d "_" `

echo $sample


samtools collate -o "$sample"_namecollate.bam "$sample"_sort.bam
done

for file in *_namecollate.bam
do
sample=`echo $file | cut -f1,2,3,4,5 -d "_" `
echo $sample

samtools fixmate -m "$sample"_namecollate.bam ../5_filtered_bam_files/alt_samtools/"$sample"_sortfixmate.bam -O bam

done

cd ../5_filtered_bam_files/alt_samtools/

for file in *_sortfixmate.bam
do
sample=`echo $file | cut -f1,2,3,4,5 -d "_" `
echo $sample

samtools sort "$sample"_sortfixmate.bam -o "$sample"_fixmate.bam


done

for file in *_fixmate.bam
do
sample=`echo $file | cut -f1,2,3,4,5 -d "_" `

samtools markdup -r "$sample"_fixmate.bam "$sample"_sortfltr.bam


done


mkdir ./sort_filt_flagstat
for file in *_sortfltr.bam

do

sample=`echo $file | cut -f1 -d "." `

echo $sample
samtools flagstat "$sample".bam \
> ./sort_filt_flagstat/"$sample"_flagstat.txt \
2> ./sort_filt_flagstat/"$sample"_flgstbam.stderr

done


cd ./sort_filt_flagstat/
#this is just combining all the bam stats files for each sample into one text file that is tab delimited so one sample per row so can then easily manipulate to summarize mapping stats etc

for filename in *_sortfltr_flagstat.txt; do
    cat "$filename"|tr '\n' '\t'
    echo "$filename"
done > ./All_RAD_sort_filt_combined.flagstat_tabbed.txt

#few lines to generate a cleaned up summary file of mapping stats post filtering

awk '{print NF}' All_RAD_sort_filt_combined.flagstat_tabbed.txt #88 columns
awk '{print $88 "\t" $1 "\t" $23 "\t" $27 "\t"$30 "\t"$44 "\t" $49 "\t" $52 "\t"$60 "\t" $64 "\t" $67 "\t" $77 }' All_RAD_sort_filt_combined.flagstat_tabbed.txt > All_RAD_sort_filt_combined_flagstat_reformat.txt
sed -i 's/(//g' All_RAD_sort_filt_combined_flagstat_reformat.txt #-i changes it in the original file; could test sending to another outfile first
sed -i 's/%//g' All_RAD_sort_filt_combined_flagstat_reformat.txt
cat ../../../scripts_and_keyfiles/short_flagstat_headers All_RAD_sort_filt_combined_flagstat_reformat.txt >All_head_filt_flagstat_reformat.txt
