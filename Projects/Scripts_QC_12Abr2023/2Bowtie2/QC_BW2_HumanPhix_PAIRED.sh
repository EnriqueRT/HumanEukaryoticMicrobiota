#!/bin/bash

#1. Get Initial Parameters
INPUT_DIR=$1
OUTPUT_DIR=$2
BW2_DB=$3

#2. Echo paths
echo "Bowtie2 QC PE:"
echo "INPUT_DIR:" $INPUT_DIR
echo "OUTPUT_DIR:" $OUTPUT_DIR
echo "BW2_DB:" $BW2_DB

#3. Filter Host-Contaminant Reads using Bowtie2 + Samtools (PAIRED FILES)
##Loop for processing samples
for f1 in $(ls $INPUT_DIR/*_1.fastq.gz)
do
	##Set name and path to pair FASTQ REVERSE
	base=$(basename $f1 _1.fastq.gz)
	f2=${INPUT_DIR}/${base}_2.fastq.gz

	##Output names
	sam_file=${OUTPUT_DIR}/${base}".sam"
	bam_file=${OUTPUT_DIR}/${base}".bam"
	bam_ffile=${OUTPUT_DIR}/${base}"_unmapped.bam"
	bam_ffile_sorted=${OUTPUT_DIR}/${base}"_unmapped_sorted.bam"
	f1_qc=${OUTPUT_DIR}/${base}"_fHP_1.fastq.gz"
	f2_qc=${OUTPUT_DIR}/${base}"_fHP_2.fastq.gz"

	##Show message
	echo "Processing Sample..."
	echo $base

	##Map with Bowtie2
	## - 8 / run in parallel with 8 threads
	echo "Mapping against Reference(Bowtie2)..."
	bowtie2 -p 8 -x $BW2_DB -1 $f1 -2 $f2 -S $sam_file

	##Convert sam to bam file
	## -@ 8 / run in parallel with 8 threads
	echo "Converting to bam file..."
	samtools view -@ 8 -bS $sam_file > $bam_file
	##Remove Sam file
	echo "Removing sam file..."
	rm $sam_file

	##SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
	## -@ 8 / run in parallel with 8 threads
	##-f 12     Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
	##-F 256   Do not(-F) extract alignments which are: <not primary alignment>
	echo "Get unmapped pairs..."
	samtools view -@ 8 -b -f 12 -F 256 $bam_file > $bam_ffile
	##Remove bam file
	echo "Removing bam file..."
	rm $bam_file

	##sort paired read alignment .bam file (sort by name -n)
	## -@ 8 / run in parallel with 8 threads
	echo "Sorting bam filtered file..."
	samtools sort -n -@ 8 $bam_ffile -o $bam_ffile_sorted
	##Remove bam unmapped-unsorted file
	echo "Removing filtered unsorted bam file..."
	rm $bam_ffile

	##split paired-end reads into separated fastq files
	## -@ 8 / run in parallel with 8 threads
	## -n / Using -n causes read names to be left as they are. 
	## -0 /dev/null | -s /dev/null | Discarding singletons, supplementary and secondary reads.
	echo "Create Separated FASTQ files..."
	samtools fastq -@ 8 $bam_ffile_sorted \
    	-1 $f1_qc \
    	-2 $f2_qc \
    	-0 /dev/null -s /dev/null -n
	##Remove bam unmapped-sorted file
	echo "Removing filtered sorted bam file..."
	rm $bam_ffile_sorted
done

#End script
echo "End Script(QC_BW2_HumanPhix_PAIRED.sh)!"
