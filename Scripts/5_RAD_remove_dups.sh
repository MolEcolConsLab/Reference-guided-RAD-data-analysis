#!/bin/bash
#Remove PCR duplicates from sorted bam files, then generate flagstats
#run from within directory containing sorted bam files
#BSUB -q short
#BSUB -W 4:00
#BSUB -R rusage[mem=1600]
#BSUB -n 6
#BSUB -R span[hosts=1]
#BSUB -e RADe_rmdups.err
#BSUB -oo RAD_rmdups.log


module load samtools/1.9
mkdir ../5_filtered_bam_files
for file in *_sort.bam

do

sample=`echo $file | cut -f1,2,3,4 -d "_" `

echo $sample

samtools rmdup "$sample"_sort.bam ../5_filtered_bam_files/"$sample"_sortfltr.bam

done

wait

cd ../5_filtered_bam_files/
mkdir ./sort_filt_flagstat
for file in *.bam

do

sample=`echo $file | cut -f1 -d "." `

echo $sample
samtools flagstat "$sample"_sortfltr.bam \
> ./sort_filt_flagstat/"$sample"_sortfltr_flagstat.txt \
2> ./sort_filt_flagstat/"$sample"_sortfltr_flgstbam.stderr

done

wait

cd ./sort_filt_flagstat/
#this is just combining all the bam stats files for each sample into one text file that is tab delimited so one sample per row so can then easily manipulate to summarize mapping stats etc

for filename in *_sortfltr_flagstat.txt; do
    cat "$filename"|tr '\n' '\t'
    echo "$filename"
done > ./All_RAD_sort_filt_combined.flagstat_tabbed.txt

wait
#few lines to generate a cleaned up summary file of mapping stats post filtering

awk '{print NF}' All_RAD_sort_filt_combined.flagstat_tabbed.txt #88 columns
awk '{print $88 "\t" $1 "\t" $23 "\t" $27 "\t"$30 "\t"$44 "\t" $49 "\t" $52 "\t"$60 "\t" $64 "\t" $67 "\t" $77 }' All_RAD_sort_filt_combined.flagstat_tabbed.txt >All_Rapture_flagstat_reformat.txt
sed -i 's/(//g' All_RAD_sort_filt_combined_flagstat_reformat.txt #-i changes it in the original file; could test sending to another outfile first
sed -i 's/%//g' All_RAD_sort_filt_combined_flagstat_reformat.txt
cat ../short_flagstat_headers All_RAD_sort_filt_combined_flagstat_reformat.txt >All_head_filt_flagstat_reformat.txt

done
