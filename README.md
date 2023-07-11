# Human Gut Eukaryotic Microbiota Project

This repository contains all the programs and scripts used to generate the results of this project. The aim of this study is to explore a common group of eukaryotic taxa (core) in the microbiota of healthy individuals across different human populations worldwide. To this end, a collection of 33 public studies and meta-studies has been compiled, including the Human Microbiome Project (HMP) and Metagenomics of the Human Intestinal Tract (MetaHIT). 

## Software

- Omics Dataset Curation Toolkit (v.1.0)
- Fastp (v.0.23.2)
- Bowtie2 (v.2.5.1)
- Samtools (v.1.16.1)
- FastQC (v.0.12.1)
- MultiQC (v.1.14)
- EukDetect (v.1.3)
- Snakemake (v.5.3.0)
- R (v.4.2.2)
- Python3 (>= 3.10)

## Python libraries

- matplotlib
- argparse
- fnmatch
- seaborn
- pandas
- numpy

## R libraries

- curatedMetagenomicData
- ComplexHeatmap
- RColorBrewer
- gridExtra
- tidyverse
- reshape2
- ggplot2
- dplyr
- readr
- DT

## Repository content

The `GLOBAL_METADATA.tsv` file contains the global metadata table for all projects. The set of variables chosen are: `file_sample_name`, `run_accession`, `sample_accession`, `study_accession`, `study_code`, `SampleID`, `instrument_model`, `DNA_extraction_kit`, `individual`, `country`, `region`, `location`, `island`, `continent`, `country_HDI_group`, `host_country_HDI_group`, `host_country_HDI_group`, `host_country_HDI_group`, `host_country_HDI_group`, `continent`, `country_HDI`, `country_HDI_group`, `host_lifestyle`, `host_sub-lifestyle`, `host_state`, `host_age`, `host_age_group`, `host_sex`, `host_BMI`, `host_diet`, `host_ethnicity`, `familyID`, `family_level`, `sibling_twins`, `twin_type` and `pregnancy_status`.

### Correlation-Prevalence_Study

In this directory can be found the programs necessary to generate the final analyses of the work:

- The `matrix_correlation.ipynb` script has been used to generate the genus-genus correlation matrix.
- The `eukdetec_prevalence_genus_study.Rmd` script has been used to generate the genus level prevalence heatmaps
- The `eukdetec_prevalence_species_heatmap_study.Rmd` script has been used to generate the species level prevalence heatmaps
- The `eukdetec_prevalence_species_stackedplot_study.Rmd` script has been used to generate three stacked bar charts with the individual level prevalences

### EukDetect

In order to be able to run EukDetect, it is necessary to generate a configfile previously, customised for each project in the study. In this way, the program `generate_configfile.py` builds the necessary configfile for each case. This program needs the following input arguments to generate the file:

| Parameter | Description | 
|   :---    |    :---     |
| `-h, --help` | Show help message and exit. |
| `-i, --fastq_directory` | Bowtie Samples Data Path. Indicate rute to processed Fastqs directory. | 
| `-l, --library_layout` | Library Layout. Indicate whether reads are paired (true) or single (false). Permitted options are {true, false}.|
| `-r, --read_length` | Read Length. Indicate length of the reads. |
| `-o1, --output_dir_eukdetect` | Output Directory EukDetect Analysis. Indicate the path to the output directory for EukDetect analyses. |
| `-o2, --output_dir_configfile` | Output Directory Configfile (Optional). Indicate the path to the output directory of configfle. Output file will be created in the current directory if not indicated. |

The file `configfile_example.txt` contains an example output of this program.

The `pres-abs_matrix_and_taxonomy.py` program performs an output analysis of EukDetect and generates presence/absence table of all taxa identified for each sample and a taxonomy table with each TaxID associated with its respective rank. This program needs the following input arguments to generate the outputs:

| Parameter | Description | 
|   :---    |    :---     |
| `-h, --help` | Show help message and exit. |
| `-i, --root_path` | Root Path. Root directory from which you want to start looking for all EukDetect outputs. Make sure the EukDetect outputs are inside that directory. |
| `-e, --euk_dir_name` | EukDetect Ouput Directory Name. Indicate the name of directory where the EukDetect outputs are saved. |
| `-r, --rank` | Taxonomic Rank. Specify one of the two allowed taxonomic ranges for which you want to obtain results. Permitted options are {genus, species}.|
| `-p, --output_name_prefix` | Output Name Prefix (Optional). Indicate prefix name for the output files. |
| `-o, --output_dir` | Output Directory (Optional). Indicate the path to the output directory. Output files will be created in the current directory if not indicated. |

The outputs generated are the following files:

- Presence/absence tables at genus and species level: `presence_absence_genus_table.tsv` and `presence_absence_species_table.tsv`.
- Taxonomy tables at genus and species level: `taxonomy_genus_table.tsv` and `taxonomy_species_table.tsv`.

The `EukDetect_example.batch` script shows an example of the complete execution of the above programs in a `.batch` file which is executed in Garnatxa cluster as follows:

```
sbatch EukDetect.batch
```

### Projects

The most relevant information on the curating and processing of all 33 studies has been included in this folder. Each project contains the same information:  

- The `multiqc_stats` directory contains three tables with the statistics of the reads:

   *  `counts_track_table.tsv` is a table where the rows correspond to the samples and the columns to the reads for each sample of the three QC steps. 

   *  `percentages_track_table.tsv` is a table with the same structure, adding for each sample the percentages of reads left over with respect to the total in each of the steps with fastp and bowtie. 

   *  `qc_stats_table.xls` is a table that shows overall statistics of reads per project, the length of the shortest and longest reads, the total number of reads, the average number of reads per sample and percentages associated with the average loss of reads per sample.

- The `Metadata_Final.xls` table contains the individual metadata.

- The three HTML documents: `multiqc_report_bowtie.html`, `multiqc_report_fastp.html` and `multiqc_report_raw.html` contain the statistics provided by MultiQC for each QC step (raw, Fastp and Bowtie2), respectively.

- The `samples_description.odt` document contains a brief description of the sample curing process.

Apart from these folders, two scripts have been included in this directory:

- The `curatedMetagenomicData_example.rmd` script is an example of how the metadata tables have been obtained from the R package curatedMetagenomicData.

- The `frequency_plots.Rmd` script has been used to obtain distribution plots of the percentage of individuals/samples by continent, HDI level, age group, lifestyle, gender and health status.

### Quality_Control

In this directory can be found three folders belonging to each of the quality control steps:

- The `Fastp` folder contains the required scripts to perform the cleaning of the reads:
   *  The `fastp_qc_PAIRED.sh` and `fastp_qc_SINGLE.sh` scripts are used to clean up the `paired-end` and `single-end` file reads, respectively.
   *  The scripts `fastp_qc_PAIRED_polyg.sh` and `fastp_qc_PAIRED_polyg.sh` are used to clean and trim reads that have poly-G tails at the 3' end (when the data came from NovaSeq or NextSeq sequencers), for both `paired-end` and `single-end` files.
   *  The `fastp_quality_control_example.batch` script shows an example of the complete execution of the above scripts in a `.batch` file which is executed in Garnatxa cluster as follows

```
sbatch fastp_quality_control.batch
```

- The `Bowtie2` folder contains the required scripts to remove the DNA belonging to the host (Human Genome GRCh38.p14) and the added Illumina contaminant (Phix Phage Genome GCF_000819615.1):
   *  The `host_contamination.Rmd` script is used to download these genomes and index them with Bowtie2.
   *  The `QC_BW2_HumanPhix_PAIRED.sh` and `QC_BW2_HumanPhix_SINGLE.sh` scripts are used to give the execution commands to Bowtie2 and generate the FASTQ files from the BAM files sorted with Samtools.
   *  The `bw2_human_phix_decontamination_example.batch` script shows an example of the complete execution of the above scripts in a `.batch` file which is executed in Garnatxa cluster as follows:

```
sbatch bw2_human_phix_decontamination.batch
```

- The `Check` folder contains the required scripts to evaluate the quality of the files with FastQC and MultiQC:
   *  The `check_qc_programs_raw_example.batch`, `check_qc_programs_fastp_example.batch` and `check_qc_programs_bowtie_example.batch` scripts show an example of the execution of FastQC and MultiQC for each of the QC steps (raw, Fastp and Bowtie2) in a `.batch` file which is executed in Garnatxa cluster as follows:

```
sbatch check_qc_programs.batch
```

### Sequencing_Depth_Filter

- The `Final_Track_Tables` folder contains all `counts_track_table.tsv` tables for all projects and the concatenated table of all of them, `concat_track_table.tsv`, which combines all samples and their associated reads in the three QC steps.

- The `sequencing_depth_plot.Rmd` script has been used to set the sequencing depth threshold.
