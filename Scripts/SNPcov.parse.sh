#!/bin/bash

for file in *.bed

do

sample=`echo $file | cut -f1 -d "." `

echo $sample


awk -v sample2=`basename $file| cut -f1 -d "." ` '{print $1,$2, $3, sample2}' "$sample".pos20.bed > ../"$sample".pos20.bed

done
