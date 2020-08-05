#!/bin/bash -l

list=$1
#enter metadata file as argument 1 in command line (BestRADBarcodesCorrected_GTHI_metadatakey.txt)
x=1
while [ $x -le 96 ]
do

	string="sed -n ${x}p ${list}"
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1, $2, $3, $4, $5, $6}')
	set -- $var
	c1=$1 #well barcode
	c2=$2 #RAD well ID
	c3=$3 #Location
	c4=$4 #Sex
	c5=$5 #sampleID
	c6=$6 #extraction_type

#run with just RA and RB files in different subfolders; run for RA, then edit accordingly and run for RB
#run with echo to check is correct, then if all looks good remove and run without
echo mv RA/RAD-GTHI_RA_GG${c1}TGCAGG.fastq RA_renamed_files/RAD-GTHI_${c2}_${c3}_${c4}_RA.fastq

#RAD-GTHI_RA_GGAACCGAGATGCAGG.fastq wA01_Hawaii_U_RA.fastq
	x=$(( $x + 1 ))

done
