# Reference-guided RAD analysis workflow

This workflow should be used for RAD-sequencing data where there is a reference genome available for the species of interest. **The goal of this workflow is to generate a set of bam (binary sequence alignment) files from your sequencing data that are clean (only true sequences remain) and aligned to the same reference genome.** These bam files can then be used for downstream analyses, such as genotyping, SNP discovery, or population structure analysis, among other applications (these analyses are not covered in this workflow.) This workflow assumes that you are starting analysis with the raw fastq files from the sequencer. Additionally, we assume that libraries are generated using the BestRAD library prep following MEC lab adapted protocol from Ali et al., though could be modified for other RAD library preps (adjusting barcodes, dup removal steps, trimming steps etc). If a reference is unavailable, use the [Stacks de novo workflow.](https://github.com/MolEcolConsLab/Stacks)

**A note on reference genomes and RAD analysis:**
Having a reference genome allows you to answer additional questions and conduct more analyses than could be done without a reference genome, especially when your reference is annotated. [This recent paper discusses the utilty of reference genomes to conservation questions. ](https://www-ncbi-nlm-nih-gov.silk.library.umass.edu/pmc/articles/PMC6895880/). There are an increasing number of high quality reference genomes for a wide variety of taxa, such as those created as part of the [Vertebrate Genomes Project](https://vertebrategenomesproject.org/). RAD-seq (restriction enzyme associated sequencing) is often chosen for genomics work in non-model species because it is a cost-effective way to answer population genomics questions and doesn't require extensive genomic resources, but [there are also caveats and limitations to what we can do with RAD-seq data.](https://www.molecularecologist.com/2017/04/17/to-radseq-or-not-to-radseq/)

**A note on scripts in this workflow:**
All scripts (bash and R scripts) are located in this [project repository](https://github.com/MolEcolConsLab/Reference-guided-RAD-data-analysis/tree/master/Scripts). To use these you should have access to and working knowledge of our cluster (for more information see [intro to MGHPCC](https://gist.github.com/MolEcolConsLab/540e2331bfe1a1147ebf0d6fc875e513) and [lab guidelines](https://hackmd.io/DkKnjNHAQmyGI98exzyFlg).)They are numbered in order of use in the workflow. Unless otherwise indicated below, all scripts contain job submission parameters for our cluster (the #BSUB commands at top of scripts). These can be modified based on cluster performance data in the future, but they are what we have used for previous projects. Projects with more/less data will require adjusting these parameters accordingly.

Occasionaly in this workflow I will make a comment about job submission or items specific to our cluster. These notes will be bulleted and begin with "C.N." 

To submit scripts as jobs, enter the following from the command line:
`[User@ghpcc06 ~]$ bsub <./script/to/run `
You can check on jobs using the `bjobs` command



**Broad workflow overview:**
![](https://i.imgur.com/d1o7gt9.png)


## 1. Metadata and raw data
### File structure:
Create one directory for the project, named with species latin name abbreviation (e.g. "chmy" for chelonia mydas) and a project title.

**To begin, create the following directories:**
1. scripts_and_keyfiles
1. 1_raw_fastq
1. 2_demultiplexed_renamed_fastq
1. 3_raw_data_quality
    1. (optional)trimmed_fastq
    1. (optional)trimmed_data_quality
1. 4_mapped_bam_files
1. 5_filtered_bam_files

**Files you need to begin:**
1. raw reads in fastq file format- these need to be unzipped prior to beginning workflow
2. metadata key file previously created during library prep linking RAD barcodes with sample IDs. This should include:
	* 	Sample ID
	* 	well ID (e.g. wA01,wB01)
	* 	Sample identifying information (e.g. species, sample collection date, sample type)
	* 	**RAD barcodes for each well** 
		* 	This is **CRUCIAL** as samples will be demultiplexed by barcode. 
		* 	The barcode list is located at: Box\MEC_lab_shared_resources\Lab protocols\RAD-seq&GBS\BestRADBarcodesCorrected_positionkey.txt
		* 	**DO NOT MODIFY THE BARCODE TEXT FILE.** Be careful with sorting spreadsheets to ensure that barcode sequences match your plate wells. Double check this step. 

* An example key file for demultiplexing can be found in this github project repository.
    * It is important that you know what column headers are for your metadata, but the column headers need to be removed for demultiplexing script to run
    * note that if you have a blank well that was sequenced, fill in metadata columns with 'NA' without any backslash or hash, this will impede the renaming script downstream

**Additional steps before beginning analyses:**
1. Clone this github repository containing all scripts to your local computer and upload them to your cluster space (can upload with filezilla, or directly clone the repository to your cluster space)
    * If you see that scripts or keyfiles have ^M in any of them, you might need to change your line endings- open these in a text editor like Atom to make sure you have unix line endings or your scripts won't run
3.  Add keyfile to the scripts_and_keyfiles directory 
4. create hackMD, github repo, or other way to document steps taken in analyses
5. Each project should also have a seperate readme doc, such as a hackmd for recording project information and workflow. 
6. Check the versions of software called in the required scripts. Make sure these are up to date, and if updates are required email the cluster administrators and update scripts as required. 
7. Check that there is enough space available in our shared project space before you start analysis. 
8. Make sure that the md5sum matches between the ftp download for the file and the files on the project space in the cluster.
	* from command line, run: `[user@ghpcc06 1_raw_fastq]$ md5sum *.fastq >raw_fastq_md5sums.txt` This will run the md5sum check for all fastq files in your current directory and output them to a text file. Each row will have the file name and unique md5sum.


## 2. Demultiplex and rename raw sequence files

**Purpose:**
Parse out pooled sequences into sample-specific sequences. If you included multiple plates on the same sequencing run, you may need to first demultiplex by plate id, and then by well id. This step will also flip and trim reads. Reads with duplicate or missing barcodes will be dropped. During library prep we use a unique barcode to identify each sample in each well, so here we will parse out sequences by each barcode and rename files with metadata from our key file. 

**Associated scripts:** 
* *BarcodeSplitListBestRadPairedEnd.pl*: this script is called by the 2a_run_BestRadSplit.sh script to parse sequences by unique barcode, flip, and trim sequences
* *2a_runBestRadSplit.sh* : This script contains fastq file locations, list of unique barcodes, and outputs demultiplexed files
* *2b_RADrename.sh* (optional) : This script renames demultiplexed barcodes with sample identifier information, but can sometimes be better to leave samples with "blind" identifiers to limit bias in analysis.
	* This script will need to be modified to create informative sample names specific to your project and number of files you have




**i. Demultiplex data by plate** 
The illumina barcodes that we add at the last library prep step (using NEB Ultra DNA library prep kit) are the plate barcodes in this protocol. So this step may be done by the sequencing center (Novogene does this for example, UCD Genome core doesn't), and is only necessary if you sequenced more than one plate a time. If you need to do this for your project, ask Lisa for the script and more details.

**ii. Demultiplex data by barcode.** 
* The 2a_runBestRadSplit.sh script should be run from within the raw_fastq directory. 
* In the script, change the following items in accordance with your project file locations:

```
f1=/project/uma_lisa_komoroske/PROJECT_DIRECTORY/raw_fastq/RAW_DATA_FILE_NAME_1.fq
f2=/project/uma_lisa_komoroske/PROJECT_DIRECTORY/raw_fastq/RAW_DATA_FILE_NAME_2.fq
out=RAD-PROJECTNAME
```
* Make sure to also change the path to the BarcodeSplitListBEstRadPairedEnd.pl script if it is not in your scripts_and_keyfiles directory. Specify out= to a descriptive project name, like RAD-GTHI (incorporates species and location)


**Expected output:**
2 Files per sample, labeled with RA and RB and the sample barcode, like so:
`RAD-GTHI_RA_GGAACCGAGATGCAGG.fastq ;
RAD-GTHI_RB_GGAACCGAGATGCAGG.fastq`



**iii. Rename files using the metadata key** (optional but recommended)
* Create two subdirectories in the 1_raw_fastq directory, one for all RA files and another for RB files,  Move files appropriately into RA and RB files
* Modify the 2b_RAD_rename.sh script; 
	* Modify each field in the following section to match your metadata key (NOTE that this is an example of some metadata that we might include in our sample names, and should be adjusted for your project)
	* You can only use the first 9 columns of your metadata key, if you have important metadata in other columns you should rearrange so that they are in the first 9 columns or use different syntax for variables.
	
	```
		var=$(echo $str | awk -F"\t" '{print $1, $2, $3, $4, $5, $6}')
			set -- $var
			c1=$1 #well barcode
			c2=$2 #RAD well ID
			c3=$3 #Location
			c4=$4 #Sex
			c5=$5 #sampleID
			c6=$6 #extraction_type

	``` 






	* If your barcode sequences are not in column 1 of your metadata key, change above to the correct column, same goes for column 2 with the well ID 
	* Variables c3-c6 are optional and can be changed based on your project metadata
	* Modify the part of the script that says `while [ $x -le 96 ]` and change 96 to the total number of samples that you have for the project


* First run the script with the last line of code hashed out to make sure your naming scheme is working properly before actually using the mv command, you should see in your log file that the names you want are actually showing up:
```
#RAD-GTHI_RA_GGAACCGAGATGCAGG.fastq wA01_Hawaii_U_RA.fastq
        x=$(( $x + 1 ))
```


* *C.N.* To run the script, don't submit as a job because this runs very quickly, and bsub doesn't always play nice when using a script and an argument following the script. Your metadata key file should be the argument following the script. First change script permissions by running from the command line `chmod u+x 2b_RAD_rename.sh` Run the script from within the scripts_and_keyfiles directory like so:
`[user@ghpcc06 sscripts_and_keyfiles]$ ./2b_RAD_rename.sh metadatakey.txt`

* Move all RA and RB files into one folder (2_demultiplexed_renamed_fastq). Leave the files unzipped for now.

**Expected output:**
Individual RA and RB files for each sample with unique identifiers, for example `MHC_Sample25_lk_MA_RA.fastq ; MHC_Sample25_lk_MA_RB.fastq`


## 3. Quality checks of raw data
**Purpose:** check initial data quality to identify failed samples, make sure files aren't corrupted, and identify adaptor contamination. 

**Associated scripts:**
* *3a_fastqcRA.sh* : Runs fastqc program on all forward reads
* *3a_fastqcRB.sh* : Runs fastqc program on all reverse reads, these two scripts can be submitted as jobs simultaneously to increase efficiency. 
* *EDA.R script, section I* : Gather additional summary stats of sample coverage
* MAKE SURE that fastq files are unzipped at this point


**i. Check line counts and file sizes**
* RA and RB files for the same sample should have the same number of lines, and should have the same file size.
* Get line counts and output to a text file. This will ONLY work on unzipped fastq files. From the directory that contains your demultiplexed fastq files (2_demultiplexed_renamed_fastq) run from the command line run:
```
du -h *fastq | sort -h >filesize_bysample.txt

wc -l *fastq >linecount_bysample.txt
```
* After getting these two text files, gzip files to save space. Run a job from the command line to do this:

`bsub -q long -W 8:00 -R rusage[mem=16000] -n 4 -R span\[hosts=1\] gzip *.fastq` 

* You should also calculate reads lost due to missing barcodes at this step. Do this by looking at the Novogene (or other sequencer) report that indcates the total number of reads generated. Then look at your linecount_by_sample.txt document, there is a total line count at the bottom. Divide this value by 4 (each sequence has 4 lines on fastq files) to get the total number of reads retained after demultiplexing. We usually see this is around 85% retained, or 15% of reads lost due to missing barcodes


**ii. Run fastqc scripts**
First run the 3a_fastqc_RA and 3a_fastqc_RB.sh scripts from the demultiplexed_renamed_fastq directory.
*  If your files are zipped make sure to change these scripts to say `*.fastq.gz` instead of `*.fastq`
```
bsub <./3a_fastqc_RA.sh
bsub <./3a_fastqc_RB.sh
```

* **Expected output:**
4 files per sample, for example: 
``
RAD-GTHI_wH11_Oahu_0_RA_fastqc.html ;
RAD-GTHI_wH11_Oahu_0_RA_fastqc.zip ;  
RAD-GTHI_wH11_Oahu_0_RB_fastqc.html ;
RAD-GTHI_wH11_Oahu_0_RB_fastqc.zip
``


	* You can look at the html files individually to see stats for each sample. Documentation for understanding all stats is [here](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)
* Move all files into the 3_raw_data_quality directory. 



**iii. Run multiqc**
* First load the multiqc module, then in the following line run multiqc from within the 3_raw_data_qualtiy directory
	* *C.N.* Note that this runs very quickly and it is not necessary to submit as a job
	* The --interactive flag allows for interactive plotting in html, so not required but nice to be able to highlight samples by metadata in the plots.

```
module load python3/3.5.0_packages/multiqc/1.4
multiqc ./ --interactive
#the multiqc command runs the program on all fastqc files in your current directory
```



* **Expected output:**
After running multiqc you should see a subdirectory called "multiqc_data" and a file called "multiqc_report.html"

* Use filezilla or another FTP client to download your multiqc report and open it in a web browser. You can also see some general stats in the multiqc_data folder when looking at "multiqc_general_stats.txt". 
* **Important things to note in the multiqc report:**
	* Only some of the metrics are really useful for our specific type of data and its important to know what to expect in RAD (vs. RNA-Seq or WG shotgun data etc)
    * Sequence duplication levels: *if very high, this may indicate high PCR duplication which will be removed in step 5 of this workflow, where you may filter out many sequences during that step. HOWEVER, in RAD data we are trying to get higher coverage of the same RAD-tag loci (remember we are pulling down a small % of the genome) that all start at the same exact cut site position, and this can results in more 'duplication' according to FastQC than other sequencing methods so don't panic if this is somewhat high...we'll really see what is going on after Step 5*
	* Mean quality scores: *These should be between 30-40 to indicate high quality. If quality scores fall below the green threshold, sequences can be trimmed for quality (you can also filter for quality downstream)*
	* Adapter contamination: *If adaptor contamination is high, sequences should be trimmed to remove adaptor sequences. Having small leftover % is generally ok and doesn't interfere with mapping*
	* For more information and explanation of other stats, see the above fastqc documentation


**iv. Generate additional summary statistics in R**
* Use the EDA.R script, section I "Initial QC check of raw read data"
* Download the filesize_bysample.txt and linecount_bysample.txt files and pull into R
* With this script, get mean, min, max, and median  number of reads by plate (optional) and well.
* produce histograms of raw reads by plate and well 
* Generate barplots of read coverage grouped by sample metadata
	* for example, can see how number of sequences varies by sex, location, collection date






## 4. Map to reference genome and produce BAM files

**Purpose:**
Align reads (fastq/fastq.gz files) that passed initial qc tests to a reference genome. Following mapping, sort sam files so that all sequences are in the same order as the genome, and convert output sam files to bam (binary) file format. The bwa (burrows wheeler aligner) is a widely used alignment algorithm. We will use the newest and fastest algorithm version, bwa-mem. [More information on this aligner here.](http://bio-bwa.sourceforge.net/)

We will use samtools to sort sam files after mapping and convert to bam files. [Samtools documentation is here.](http://www.htslib.org/doc/samtools.html) Samtools sorts sam files based on leftmost position of alignments to the reference genome. After generating bam files, the script calls  a function in samtools again to get summary statistics of mapping for each sample.

**Associated scripts:**
*4_RAD_bwa.sh* : indexes reference genome, maps reads, sorts SAM files, converts sam to bam, and generates individual mapping stats


**i. Prepare reference genome**
* create a directory called "reference_genome"n your main project directory
* Download the genome in fasta format from NCBI and add to this directory
	* Make sure to indicate which reference was used and the version as well as access date in your project readme file. 
	* Find your reference genome on the NCBI website, for example I've used this one, the file should end in .fna.gz: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Arvicola_amphibius/all_assembly_versions/GCA_903992535.1_mArvAmp1.1/GCA_903992535.1_mArvAmp1.1_genomic.fna.gz
	* To download to the cluster, navigate to your reference_genome directory, and from the command line run: `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Arvicola_amphibius/all_assembly_versions/GCA_903992535.1_mArvAmp1.1/GCA_903992535.1_mArvAmp1.1_genomic.fna.gz`
	    * replace this link to your appropriate reference genome


**ii. Align reads to reference genome**
* You will run the 4_RAD_bwa.sh script from the directory that contains your qc-passed fastq files.
	*  Make sure you have the "short_flagstat_headers" file in your scripts_andkey_files directory prior to running this script. 
* Script components to modify: 
	* the `bwa index` command indexes the reference genome, which is needed prior to alignment. Change the path to point to the location of your reference genome
	* If you gzipped your fastq files above, change `for file in ./*_RA.fastq` to `for file in ./*_RA.fastq.gz`
	* Change `sample=echo $file |cut -f1,2,3,4 -d "_"`
 line of the script based on how many fields are in your sample names. (e.g. the file "RAD-GTHI_wC08_EastFFS_M_RB.fastq" has 5 fields delimited by _ and I want my sample name to contain the first 4 fields, before the RB)
	* *Flags used in script:*
		* 	-t : this is the number of threads to use for the job, which should match the -n flag in #BSUB commands
		*   -O : from samtools sort, this indicates the output file should be in bam formatis output in bam format
		* -o : specifies the output file name. This should be "sample"_sort.bam. We also generate standard error files for each sample. 
* *What the script is doing:*
	* maps all reads to the reference genome with `bwa mem` command
	* gets mapping statistics for all samples with the `samtools flagstat` command
	* Combine individual mapping statistics into a combined text file 
	
**Expected output:**
* One bam file per sample, and one stderr file per sample, ocated in the directory "4_mapped_bam_files" 
	* Since bam files are in binary format, you can't view them directly, but could use the "samtools view" command to view them.  
* One flagstat.txt file and one flgstbam.stdrr file for each sample, which will be located in the directory "sort_flagstat" within the "4_mapped_bam_files" directory. 
* Combined flagstat file called “All_head_flagstat_reformat.txt” which has  mapping stats for all samples, and reformat to include column headers and remove unnecessary columns using the short_fagstat_headers file


## 5. Remove PCR duplicates

**Purpose:**
For library preps where the final step is PCR amplification, it is common to generate PCR duplicates and sequence them, which we want to remove as including these duplicates result in biased estimates of coverage and downstream analyses. This script also generates mapping statistics for all samples post-filtering of PCR duplicates.

**Associated scripts:**
* *5_removedups.sh* : remove PCR duplicates using samtools, create a combined summary stats text file

**i. Remove PCR duplicates and generate post-removal mapping stats**
* Run the 5_remove_dups.sh script from within the directory containing sorted bam files 
	
* *What the script is doing:*
	* removes duplicates from bam files using the samtools `rmdup` command. 
	* get mapping statistics using the samtools `flagstat`
	* combine flagstats from all samples into one file, and  reformat to include column headers and remove unnecessary columns using the short_fagstat_headers file


* Alternatively, you can use the more recent version of samtools and the markdup option instead, however I found this to have the same result (removing the same amount of duplicates) as keeping the old version of samtools used in this script. 
    * *Use script 5_remove_dups_altsamtools.sh*
    * You can also use picardtools to do this step, but we do not currently have a script for it. 
	

**Expected output:**
* sorted, filtered bam files for every sample, 
*  stderr file for every sample
* Flagstat file for every sample
* Combined flagstat file called "All_head_filt_flagstat_reformat.txt" which has post-filtering mapping stats 




## 6. Examine mapping QC metrics

**Purpose:**
Compare the filtered and unfiltered mapping stats to estimate the percent loss from PCR duplicates, determine which samples fail or pass QC based on a number of metrics, and estimate coverage per locus and sample. *Results of this will help determine if good to move on to downstream steps or if need to conduct further QC or otherwise deal with issues that arise (so need to actually work through the output to understand your data, dont just generate and move on!)*

**Associated scripts:**
* EDA.R section II

**Required input files:**
* Previous raw linecount.txt file from part 3
* All_head_flagstat_reformat.txt from part 4
* All_head_filt_flagstat_reformat.txt from part 5
* File specifying DNA input concentrations per sample so you can see effects of input concentration

**What the script is doing:**

*  pull in the mapping stats files from mapping the raw sequences to your reference genome, and the mapping stats after removing PCR duplicates. 
*  get percentage of raw reads left after filtering
*  get proportion of reads that were PCR duplicates
*  use qc info from section 1 to identify failed samples
*  check the effects of DNA input concentration


## Determine coverage
We want to estimate coverage in the sequenced areas acorss samples. Lisa will add in script for this.

## 7. Index bam files (optional)
Indexed bam files are required for input into the Freebayes Genotyping workflow, viewing files in IGV, and some other workflows. This script uses samtools to index the bam files against your reference genome and allows for comparison between bam files based on their relative positions. 

**Associated scripts:**
* *index_bam.sh*

**Running the script:**
* Run the short index_bam.sh script from within the "5_filtered_bam_files" directory, which contains all of your filtered, sorted bam files
	* C.N. You don't need to specify any job submission parameters, since this runs quickly on the short queue and default job submission parameters are fine

* the $i in the script corresponds to the name of each bam file in this directory, so all bam files are going to be indexed. If you don't wan tot index all bam files, you should pull those you don't want to index into a subdirectory.

**Expected output:**
* In the "5_filtered_bam_files" you should see a corresponding file for each sample that ends in the file extension .bai.
	* Keep the .bam and .bai files in the same directory, this is necessary when you need indexed bam files in downstream workflows
* You'll get an email from the cluster where each line says "indexing: sample_name"

-------------------
### From here you can proceed to downstream workflows, such as genotyping with freebayes.