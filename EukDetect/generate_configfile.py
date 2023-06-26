#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#*****************************************************************************#
#                                                                             #
# configfile_samples.py                                                       #
#                                                                             #
# This program generates the configfile with the necessary parameters for     #
# running EukDetect pipeline using snakemake.                                 #
#                                                                             #
# Version 1.0.                                                                #
#                                                                             #
# Author: Enrique Roig Tormo                                                  #
# Date: 24/04/23                                                              #
#                                                                             #
#*****************************************************************************#

#IMPORTS
import argparse
import fnmatch
import os
import re


#FUNCTIONS
def samples_names(input_dir):
    """
    This function generates a sorted list of the unique names of all fastqs 
    files in a directory.   

    Parameters
    ----------
    input_dir : str
        The input directory where are located all fastqs files.

    Returns
    -------
    unique_names : list
        List with the unique samples names sorted.

    """

    #List files with .txt extension in the current directory
    fastqs = fnmatch.filter(os.listdir(input_dir), '*fastq.gz')

    sample_name = []

    #Obtain samples names of fastqs
    for fastq in fastqs:
        patron = r"_\d\.fastq\.gz$"
        subcadena = re.split(patron, fastq)[0]
        sample_name.append(subcadena)
    
    #Convert the list to a set to remove duplicates, then convert back to a list
    unique_name = list(set(sample_name))
    #Sort the list
    unique_name.sort()

    return unique_name


def write_configfile(fq_dir, paired_end, readlen, output_euk, output_dir, list_names):
    """
    This function writes in a new file all information that needs the 
    configfile names of the samples according to the 
    format of the config_file. If an output path has been provided, the file 
    will be saved in that path, otherwise it will be automatically saved in 
    the directory where the program is located.

    Parameters
    ----------
    fq_dir : str
    
    paired_end : str
        This variable takes the value true if the samples are paired or false 
        if they are single
    readlen : int
        Length of the reads
    output_euk : str
        Output directory 
    output_dir : str
        The provided output directory.
    list_names : list
        List with the unique samples names sorted.

    Returns
    -------
    None

    """

    #Open the new configfile in the currently directory or in a specified one.
    if output_dir:
        #Join the output path and filename to create the full path
        final_path = os.path.join(output_dir, "configfile.txt")
    else:
        current_path = os.getcwd()
        final_path = current_path + "/configfile.txt"

    with open(final_path, 'w') as f:
        f.write("#Default config file for eukdetect. Copy and edit for analysis\n")
        f.write("\n#Directory where EukDetect output should be written\n")
        f.write("output_dir: " + output_euk + "\n")
        f.write("\n#Indicate whether reads are paired (true) or single (false)\n")
        f.write("paired_end: " + paired_end + "\n")

    if paired_end == "true":
        with open(final_path, "a") as f:
            f.write("\n#filename excluding sample name. No need to edit if paired_end = false\n")
            f.write("fwd_suffix: _1.fastq.gz\n")
            f.write("\n#filename excludign sample name. no need to edit if paired_end = false\n")
            f.write("rev_suffix: _2.fastq.gz\n")
    else:
        with open(final_path, "a") as f:
            f.write("\n#file name excluding sample name. no need to edit if paired_end = true\n")
            f.write("se_suffix: .fastq.gz\n")

    with open(final_path, 'a') as f:
        f.write("\n#length of your reads. pre-trimming reads not recommended\n")
        f.write("readlen: " + str(readlen) + "\n")
        f.write("\n#full path to directory with raw fastq files\n")
        f.write("fq_dir: " + fq_dir + "\n")
        f.write("\n#full path to folder with eukdetect database files\n")
        f.write("database_dir: /home/enroig2/metagenomics/EukDetect/eukdb\n")
        f.write("\n#name of database. Default is original genomes only database name\n")
        f.write("database_prefix: ncbi_eukprot_met_arch_markers.fna\n")
        f.write("\n#full path to eukdetect installation folder\n")
        f.write("eukdetect_dir: /home/enroig2/metagenomics/EukDetect\n")
        f.write("\n#list sample names here. fastqs must correspond to {samplename}{se_suffix} for SE reads or {samplename}{fwd_suffix} and {samplename}{rev_suffix} for PE\n")
        f.write("#each sample name should be preceded by 2 spaces and followed by a colon character\n")
        f.write("samples:\n")

        for i in list_names:
            f.write("  " + i + ":\n")


#MAIN PROGRAM
def main():

    #Setting Arguments
    parser = argparse.ArgumentParser()

    ##Parameter fastq_directory
    parser.add_argument(
        '-i','--fastq_directory',
        action = 'store',
        required = True,
        help ='Bowtie Samples Data Path. Indicate rute to processed Fastqs directory.'
        )
    
    ##Parameter library_layout
    parser.add_argument(
        '-l','--library_layout',
        action = 'store',
        required = True,
        help = 'Library Layout. Indicate whether reads are paired (true) or single (false)'
        )

    ##Parameter read_length
    parser.add_argument(
        '-r','--read_length',
        action = 'store',
        required = True,
        help = 'Read Length. Indicate length of the reads.'
        )

    ##Parameter output_directory
    parser.add_argument(
        '-o1','--output_dir_eukdetect', 
        action = 'store',
        required = True,
        help = 'Output Directory EukDetect Analysis. Indicate the path to the' 
        'output directory for EukDetect analyses.'
    )

    ##Parameter output_directory
    parser.add_argument(
        '-o2','--output_dir_configfile', 
        action = 'store',
        required = False,
        help = 'Output Directory Configfile (Optional). Indicate the path to' 
        'the output directory of configfle.' 
        'Output file will be created in the current directory if not indicated.'
    )

    #Process arguments
    args = parser.parse_args()
    fastq_dir = args.fastq_directory
    lib_layout_raw = args.library_layout
    library_layout = lib_layout_raw.lower()
    read_length = args.read_length
    output_euk = args.output_dir_eukdetect
    output_dir = args.output_dir_configfile

    #Call functions
    list_names = samples_names(fastq_dir)
    write_configfile(fastq_dir, library_layout, read_length, output_euk, 
                     output_dir, list_names)


if __name__ == '__main__':
    main()