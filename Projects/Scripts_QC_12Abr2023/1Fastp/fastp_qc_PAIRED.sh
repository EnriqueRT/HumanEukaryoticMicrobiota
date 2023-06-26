#!/bin/bash

#1. Get Initial Parameters
INPUT_DIR=$1
OUTPUT_DIR=$2

#2. Echo paths
echo "Fastp QC PE:"
echo "INPUT_DIR:" $INPUT_DIR
echo "OUTPUT_DIR:" $OUTPUT_DIR

#3. Do Quality Control of reads with Fastp (PAIRED)
##Loop for processing samples
for f1 in $(ls $INPUT_DIR/*_1.fastq.gz)
do
	###Set name and path to pair FASTQ REVERSE
	base=$(basename $f1 _1.fastq.gz)
	f2=${INPUT_DIR}/${base}_2.fastq.gz

	###Output names
	f1_qc=${OUTPUT_DIR}/${base}"_qc_1.fastq.gz"
	f2_qc=${OUTPUT_DIR}/${base}"_qc_2.fastq.gz"

	###Show message
	echo "Processing Sample..."
	echo $base

	###QC with FASTP
	#Input and Output files(R1)
	#Input and Output files(R2)
	#Adapters(AUTODETECT PE)
	#Trimming parameters(LEADING:20)
	#Trimming parameters(TRAILING:20)
	#Trimming parameters(SLIDINGWINDOW:4:15)
	#Cut tail 1(-t 1) and trim reads greater than 290 bp (-b 290)
	#Default filters
	#Low complexity filter
	#Extra Filters (N=0 and lenght=60bp)
	#We do this to not save reports generated by fastp(JSON and HTML)
	#Disable Duplication Evaluation
	#Number of threads
	fastp \
	--in1 $f1 --out1 $f1_qc \
	--in2 $f2 --out2 $f2_qc \
	--detect_adapter_for_pe \
	-5 --cut_front_mean_quality 20 --cut_front_window_size 1 \
	-3 --cut_tail_mean_quality 20 --cut_tail_window_size 1 \
	--cut_right --cut_right_mean_quality 15 --cut_right_window_size 4 \
	-t 1 -b 70 \
	--qualified_quality_phred 15 --unqualified_percent_limit 40 \
	--low_complexity_filter \
	--n_base_limit 0 --length_required 60 \
	--json /dev/null --html /dev/null \
	--dont_eval_duplication \
	--thread 8
done

#End script
echo "End Script (fastp_qc_PAIRED.sh)!"
