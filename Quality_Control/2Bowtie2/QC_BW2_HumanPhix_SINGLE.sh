#!/bin/bash

#1. Get Initial Parameters
INPUT_DIR=$1
OUTPUT_DIR=$2
BW2_DB=$3

#2. Echo paths
echo "Bowtie2 QC SE:"
echo "INPUT_DIR:" $INPUT_DIR
echo "OUTPUT_DIR:" $OUTPUT_DIR
echo "BW2_DB:" $BW2_DB

#3. Filter Host-Contaminant Reads using Bowtie2 + Samtools (SINGLE FILES)
##Loop for processing samples
for f1 in $(ls $INPUT_DIR/*.fastq.gz | grep -v "_1.fastq.gz" | grep -v "_2.fastq.gz")
do
	##Set name
	base=$(basename $f1 .fastq.gz)

	##Output names
	sam_file=${OUTPUT_DIR}/${base}".sam"
	bam_file=${OUTPUT_DIR}/${base}".bam"
	bam_ffile=${OUTPUT_DIR}/${base}"_unmapped.bam"
	bam_ffile_sorted=${OUTPUT_DIR}/${base}"_unmapped_sorted.bam"
	f1_qc=${OUTPUT_DIR}/${base}"_fHP.fastq.gz"

	##Show message
	echo "Processing Sample..."
	echo $base

	##Map with Bowtie2
	## - 8 / run in parallel with 8 threads
	echo "Mapping against Reference(Bowtie2)..."
	bowtie2 -p 8 -x $BW2_DB -U $f1 -S $sam_file

	##Convert sam to bam file
	## -@ 8 / run in parallel with 8 threads
	echo "Converting to bam file..."
	samtools view -@ 8 -bS $sam_file > $bam_file
	##Remove Sam file
	echo "Removing sam file..."
	rm $sam_file

	##SAMtools SAM-flag filter: get unmapped reads
	## -@ 8 / run in parallel with 16 threads
	## -f 4 keep unmapped reads
	##-F 256   Do not(-F) extract alignments which are: <not primary alignment>
	echo "Get unmapped reads..."
	samtools view -@ 8 -b -f 4 -F 256 $bam_file > $bam_ffile
	##Remove bam file
	echo "Removing bam file..."
	rm $bam_file

	##sort read alignment .bam file (sort by name -n)
	## -@ 8 / run in parallel with 8 threads
	echo "Sorting bam filtered file..."
	samtools sort -n -@ 8 $bam_ffile -o $bam_ffile_sorted
	##Remove bam unmapped-unsorted file
	echo "Removing filtered unsorted bam file..."
	rm $bam_ffile

	##Generate fastq file
	## -@ 8 / run in parallel with 8 threads
	## -n / Using -n causes read names to be left as they are. 
	echo "Create FASTQ file..."
	samtools fastq -@ 8 $bam_ffile_sorted -n -0 $f1_qc
	##Remove bam unmapped-sorted file
	echo "Removing filtered sorted bam file..."
	rm $bam_ffile_sorted
done

#End script
echo "End Script(QC_BW2_HumanPhix_SINGLE.sh)!"
