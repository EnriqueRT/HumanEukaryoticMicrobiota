#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#*****************************************************************************#
#                                                                             #
# total_marker_cov.py                                                         #
#                                                                             #
# This program performs an output analysis of EukDetect and generates a table #
# with the total marker coverage values for each taxon associated with each   #
# sample.                                                                     #
#                                                                             #
# Version 1.0.                                                                #
#                                                                             #
# Author: Enrique Roig Tormo                                                  #
# Date: 30/05/23                                                              #
#                                                                             #
#*****************************************************************************#

#IMPORTS
import argparse
from time import time
import pandas as pd
import fnmatch
import os


#FUNCTIONS
def find_Euk_output(input_dir, euk_dir_name):
    """
    This function finds directories matching the specified EukDetect directory 
    name within the given input directory.

    Parameters
    ----------
    input_dir : str 
        Path to the input directory to search.
    euk_dir_name : str 
        Name of the EukDetect output directory to look for.

    Returns
    -------
    euk_output_dirs : list 
        A list of paths to directories that match the specified EukDetect directory name.
    """

    euk_output_dirs = []
    for root, dirs, files in os.walk(input_dir):
        if euk_dir_name in dirs:
            euk_output_dirs.append(os.path.join(root, euk_dir_name))
    
    return euk_output_dirs


def get_process_taxonomy_table(taxo_table_path):
    """
    This function reads the taxonomy table obtained with the Presence-Absence 
    program and adds an extra column with the lineage of each taxon reconstructed 
    with the way it is originally found in the EukDetect output files.

    Parameters
    ----------
    taxo_table_path : str 
        Path to the taxonomy species table generated for the Presence-Absence Table.

    Returns
    -------
    taxonomy_tab : pandas.DataFrame 
        Taxonomy table with new column of interest.
    """
    
    taxonomy_tab = pd.read_csv(taxo_table_path, sep = "\t")

    lin_total = []

    for column, row in taxonomy_tab.iterrows():
        values = []
        for col, value in row.items():
            if pd.isnull(value) or col == "TaxID":
                continue
            values.append(f"{col}-{value}")
        lin = "|".join(values)
        lin_total.append(lin)

    taxonomy_tab["lineage"] = lin_total
    
    return taxonomy_tab


def generate_dicts(dict_taxid_sam, dict_marker_sam, count_filter, count_read, count_success, count_warnings, df, run, taxo_sp_table):
    """
    This function generates the dictionary dict_taxid_sam containing all the 
    taxids associated with a sample and the dict_marker_sam containing the 
    values of the total_marker_coverage column associated with that sample.
    In addition, the files that have been successfully processed and those that 
    have not, either because they did not pass the filters required by EukDetect 
    or because they did not align with anything, are counted.

    Parameters
    ----------
    dict_taxid_sam : dict
        Dictionary to store taxid information for each sample.
    dict_marker_sam : dict 
        Dictionary to store total_marker_coverage information for each sample.
    count_filter : int 
        Counter for the number of samples failing filter requirements.
    count_read : int 
        Counter for the number of samples with empty read count files.
    count_success : int
        Counter for the number of successfully processed samples.
    count_warnings : int
        Counter for the number of successfully processed samples with some warning.
    df : DataFrame
        Dataframe containing EukDetect output file.
    run : str
        Run identifier for the current sample.
    taxo_sp_table : pandas dataframe
        The processed taxonomy species table from Absence-Presence Table

    Returns
    -------
    count_filter : int 
        Counter for the number of samples failing filter requirements.
    count_read : int 
        Counter for the number of samples with empty read count files.
    count_success : int
        Counter for the number of successfully processed samples.
    count_warnings : int
        Counter for the number of successfully processed samples with some warning.
    """

    #TaxID and lineage for taxonomy table
    taxids_ref = set(taxo_sp_table['TaxID'])
    lineages_ref = set(taxo_sp_table['lineage'])
    
    #Init keys in dictionaries for this sample
    dict_taxid_sam[run] = []
    dict_marker_sam[run] = []    

    if df.empty:
        if df.columns == "No taxa passing filter requirements.":
            count_filter += 1
        elif df.columns == "Empty read count file. Likely no aligned reads in sample.":
            count_read += 1
            
    else:
        #Go through each taxon at the species level
        for index, row in df.iterrows():
            taxid = row["Taxid"]
            marker = row["Total_marker_coverage"]
            lineage = row["Lineage"]
            
            if (lineage in lineages_ref) and (taxid in taxids_ref):
                dict_taxid_sam[run].append(taxid)
                dict_marker_sam[run].append(marker)
            elif (lineage in lineages_ref) and len(df[df["Lineage"]==lineage])==1:
                #Get taxid for linage_ref and add ref_taxid
                taxo_filt = taxo_sp_table[taxo_sp_table["lineage"]==lineage]
                ref_taxid = list(taxo_filt['TaxID'])[0]
                dict_taxid_sam[run].append(ref_taxid)
                #Add marker
                dict_marker_sam[run].append(marker)
            else:
                #Print sample name
                print(run)
                #Add warning_taxid 
                warning_taxid = 'warning_' + str(taxid)
                dict_taxid_sam[run].append(warning_taxid)
                #Add marker
                dict_marker_sam[run].append(marker)
                count_warnings += 1
        count_success += 1

    return count_filter, count_read, count_success, count_warnings


def total_cov_table(dict_taxid_sam, dict_marker_sam, unique_taxid):
    """
    This function creates an empty dataframe with the sample names as rows and 
    the TaxIDs as columns. Then, the dict_marker_sam dictionary elements are  
    traversed, combining the corresponding values of dict_taxid_sam[sample] and 
    markers in ordered pairs, thus allowing the tot_cov_sam_tab DataFrame to be 
    updated with the values of the total_marker_coverage column (without 
    percentages) of each taxon for each sample.

    Parameters
    ----------
    dict_taxid_sam : dict 
        Dictionary containing taxid information for each sample.
    dict_marker_sam : dict 
        Dictionary containing marker information for each sample.
    unique_taxid : list
        List of unique taxids

    Returns
    -------
    tot_cov_sam_tab : DataFrame
        A dataframe representing all total_marker_coverage values per sample.
    """

    tot_cov_sam_tab = pd.DataFrame(index=dict_marker_sam.keys(), columns=unique_taxid)
    for muestra, markers in dict_marker_sam.items():
        list_mark = zip(dict_taxid_sam[muestra], markers)
        for taxid, porcentaje in list_mark:
            tot_cov_sam_tab.loc[muestra, taxid] = porcentaje.rstrip("%")

    #Fill empty spaces with 0
    tot_cov_sam_tab.fillna(0, inplace=True)

    return tot_cov_sam_tab


#MAIN PROGRAM
def main():

    print("""
    ###########################################################################
    ##                                                                       ##
    ##                        total_marker_cov.py                            ##
    ##                          Date: 30/05/23                               ##
    ##                                                                       ##
    ##   Description: This program performs an output analysis of EukDetect  ##
    ##   and generates a table with the total marker coverage values for     ##
    ##   each taxon associated with each sample.                             ##
    ##                                                                       ##
    ##   Author: Enrique Roig Tormo                                          ##
    ##                                                                       ##
    ###########################################################################
    """)

    #Setting Arguments
    parser = argparse.ArgumentParser()

    ##Parameter root_path
    parser.add_argument(
        '-i','--root_path',
        action = 'store',
        type=str,
        required = True,
        help ='Root Path. Root directory from which you want to start looking for all EukDetect outputs.' 
        'Make sure the EukDetect outputs are inside that directory.'
    )

    ##Parameter euk_dir_name
    parser.add_argument(
        '-e','--euk_dir_name',
        action = 'store',
        type=str,
        required = True,
        help ='EukDetect Ouput Directory Name. Indicate the name of directory where the EukDetect outputs are saved.'
    )
    
    ##Parameter taxonomy_sp_table
    parser.add_argument(
        '-t','--taxonomy_sp_table',
        action = 'store',
        type=str,
        required = True,
        help ='Taxonomy Species Table. Indicate the path of to the Taxonomy Species Table generated for the Presence-Absence Table.'
    )
    
    ##Parameter output_file_prefix
    parser.add_argument(
        '-p','--output_name_prefix', 
        action = 'store',
        type=str,
        required = False,
        help = 'Output Name Prefix (Optional). Indicate prefix name for the output files.'
    )

    ##Parameter output_directory
    parser.add_argument(
        '-o','--output_dir', 
        action = 'store',
        type=str,
        required = False,
        help = 'Output Directory (Optional). Indicate the path to the output directory.' 
        'Output files will be created in the current directory if not indicated.'
    )
    

    #Process arguments
    args = parser.parse_args()
    input_dir = args.root_path
    euk_dir_name = args.euk_dir_name
    taxo_table_path = args.taxonomy_sp_table
    prefix_name = args.output_name_prefix
    output_dir = args.output_dir

    #Call functions

    t0 = time()

    #List of all directories with EukDetect outputs
    euk_output_dirs = find_Euk_output(input_dir, euk_dir_name)
    
    #Get taxonomy table
    taxo_sp_table = get_process_taxonomy_table(taxo_table_path)

    #Create the empty dictionaries with the TaxIDs and total_marker_coverage 
    #associated with each sample
    dict_taxid_sam = {}
    dict_marker_sam = {}
    
    #Process each directory -study
    for path in euk_output_dirs:
        #List all files in the Euk_output directory
        euk_files = fnmatch.filter(os.listdir(path), '*hits_table.txt')
        
        #Initialise the counters for the statistics of processed files
        count_filt = 0
        count_read = 0
        count_success = 0
        count_warnings = 0

        #Loop to process each file
        for file in euk_files:
            final_path = path + "/" + file
            df = pd.read_csv(final_path, sep = "\t")
            run = file.split("_")[0]
            output = generate_dicts(dict_taxid_sam, dict_marker_sam, count_filt, count_read, count_success, count_warnings, df, run, taxo_sp_table)
            count_filt, count_read, count_success, count_warnings = output[0], output[1], output[2], output[3]
            
        #Show total files and messages for the directory-study
        tot_files = count_filt + count_read + count_success
        print("\n" + path)
        print("\nTotal de archivos procesados por EukDetect:", tot_files) 
        print("Archivos que no han conseguido superar los requisitos del filtro"
              " taxonomico:", count_filt)
        print("Archivos que presentan un recuento de lecturas vacio, probablemente porque no haya "
            "lecturas alineadas en la muestra:", count_read)
        print("Archivos que se han procesado con exito:", count_success)
        print("Total Warning Taxids calls(not equivalent to total warning_taxids):", count_warnings)
    
    #Get unique list of taxids and Generate final table
    unique_taxid = list(set(sum(dict_taxid_sam.values(), [])))
    tot_cov_sam_tab = total_cov_table(dict_taxid_sam, dict_marker_sam, unique_taxid)

    #Name of the resulting tables
    if prefix_name:
        name_tab = prefix_name + "_total_marker_table.tsv"
    else:
        name_tab = "total_marker_table.tsv"

    #Save the path of result table in a specified directory or in the currently one
    if output_dir:
        path_fi = os.path.join(output_dir, name_tab)
    else:
        current_path = os.getcwd()
        path_fi = current_path + "/" + name_tab

    #Save the results
    tot_cov_sam_tab.to_csv(path_fi, sep='\t',header=True, index_label="Sample")

    tf = time()
    tot_time = tf - t0
    print("\nTiempo de ejecucion:", tot_time)
    print("\n")


if __name__ == '__main__':
    main()