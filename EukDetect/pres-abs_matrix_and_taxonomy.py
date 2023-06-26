#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#*****************************************************************************#
#                                                                             #
# pres-abs_matrix_and_taxonomy.py                                             #
#                                                                             #
# This program performs an output analysis of EukDetect and generates         #
# presence/absence table of all taxa identified for each sample and a         #
# taxonomy table with each TaxID associated with its respective rank.         #
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


def generate_dicts(dict_taxid_sam, dict_taxid_lin, count_filter, count_read, count_success, count_na_rank, df, run, rank):
    """
    This function generates the dictionary dict_taxid_sam containing all the 
    taxids associated with a sample and the dict_taxid_lin containing the all
    lineages associated that TaxID. In addition, the files that have been 
    successfully processed and those that have not, either because they did not 
    pass the filters required by EukDetect or because they did not align with 
    anything, are counted.

    Parameters
    ----------
    dict_taxid_sam : dict
        Dictionary to store taxid information for each sample.
    dict_taxid_lin : dict 
        Dictionary to store lineages information for each TaxID.
    count_filter : int 
        Counter for the number of samples failing filter requirements.
    count_read : int 
        Counter for the number of samples with empty read count files.
    count_success : int
        Counter for the number of successfully processed samples.
    count_na_rank : int:
        Counter for the number of samples with out taxa after selecting Rank.
    df : DataFrame
        Dataframe containing EukDetect output file.
    run : str
        Run identifier for the current sample.
    rank : str
        Taxonomic rank to be filtered by

    Returns
    -------
    count_filter : int 
        Counter for the number of samples failing filter requirements.
    count_read : int 
        Counter for the number of samples with empty read count files.
    count_success : int
        Counter for the number of successfully processed samples.
    """

    if df.empty:
        if df.columns == "No taxa passing filter requirements.":
            count_filter += 1
        elif df.columns == "Empty read count file. Likely no aligned reads in sample.":
            count_read += 1
        dict_taxid_sam[run] = []
    else:
        # Filter the containing rows by the provided Rank.
        df_filt = df.loc[df["Rank"] == rank]
        dict_taxid_sam[run] = []
        
        if df_filt.empty:
            count_na_rank += 1
        else:
            for index, row in df_filt.iterrows():
                taxid = row["TaxID"]
                lineage = row["Lineage"]
                dict_taxid_sam[run].append(taxid)
                if lineage not in dict_taxid_lin:
                    dict_taxid_lin[taxid] = lineage
            count_success += 1

    return count_filter, count_read, count_success, count_na_rank


def pres_abs_table(dict_taxid_sam, dict_taxid_lin):
    """
    This function creates a presence/absence table of all taxa identified for 
    each sample from the generated dictionaries. If the taxa is found for a 
    sample, a 1 is assigned, otherwise a 0.

    Parameters
    ----------
    dict_taxid_sam : dict
        Dictionary of taxid information for each sample.
    dict_taxid_lin : dict 
        Dictionary of lineages information for each TaxID.

    Returns:
    df : DataFrame
        Presence/absence table with taxIDs as rows and samples as columns.
    """
        
    #Set results df
    df=pd.DataFrame()
    
    #Get column with sample names
    df['Sample']=list(dict_taxid_sam.keys())
    
    #Get list of all taxids
    all_taxids=list(dict_taxid_lin.keys())
    
    
    #Get pres_abs col por each taxid
    for taxid in all_taxids:
        #Get temporal copy
        temp_pres_abs = dict_taxid_sam.copy()
        #For samples in copy
        for key in temp_pres_abs:
            if taxid in temp_pres_abs[key]:
                temp_pres_abs[key]=1
            else:
                temp_pres_abs[key]=0
        #Save taxid column
        df[taxid]=df['Sample'].map(temp_pres_abs)
    
    return df

  
def taxonomy_table(rank, dict_taxid_lin):
    """
    This function creates a taxonomy table with each TaxID associated with its 
    respective rank.

    Parameters
    ----------
    rank : str
        Taxonomic rank used to determine the levels in the table.
    dict_taxid_lin : dict 
        Dictionary of lineages information for each TaxID.

    Returns
    -------
    taxonomy_tab : DataFrame 
        Taxonomy table with taxID as rows and taxonomic levels as columns.
    """

    #Set the columns of the dataframe according to rank
    if rank == "genus":
        taxa = ["phylum", "class", "order", "family", "genus"]
    else:
        taxa = ["phylum", "class", "order", "family", "genus", "species"]
    
    taxonomy_tab = pd.DataFrame(index=dict_taxid_lin.keys(), columns=taxa)
    for taxid, taxon in dict_taxid_lin.items():
        #Split the string by the character "|"
        division = taxon.split('|')
        #For each taxonomic level, look up the corresponding name
        for level in taxa:
            for sub_div in division:
                if sub_div.startswith(level + '-'):
                    #Assign the corresponding name to the DataFrame
                    taxonomy_tab.at[taxid, level] = sub_div.split('-')[1]

    return taxonomy_tab


#MAIN PROGRAM
def main():

    print("""
    ###########################################################################
    ##                                                                       ##
    ##                      pres-abs_matrix_and_taxonomy.py                  ##
    ##                          Date: 30/05/23                               ##
    ##                                                                       ##
    ##   Description: This program performs an output analysis of EukDetect  ##
    ##   and generates a presence/absence table of all taxa identified for   ##
    ##   each sample and a taxonomy table with each TaxID associated with    ##
    ##   its respective rank.                                                ##
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

    ##Parameter taxonomic rank
    parser.add_argument(
        '-r','--rank', 
        action = 'store',
        type=str,
        choices=["genus", "species"],
        required = True,
        help = 'Taxonomic Rank (genus or species). Specify one of the two allowed taxonomic ranges for which you want to obtain results.'
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
    rank = args.rank
    prefix_name = args.output_name_prefix
    output_dir = args.output_dir

    #Call functions

    t0 = time()

    euk_output_dirs = find_Euk_output(input_dir, euk_dir_name)

    #Create the empty dictionaries with the samples and lineages associated with 
    # each TaxID
    dict_taxid_sam = {}
    dict_taxid_lin = {}
    
    #Process each directory -study
    for path in euk_output_dirs:
        #List all files in the Euk_output directory
        euk_files = fnmatch.filter(os.listdir(path), '*hits_eukfrac.txt')
        
        #Initialise the counters for the statistics of processed files for the directory
        count_filt = 0
        count_read = 0
        count_success = 0
        count_na_rank = 0

        #Loop to process each file
        for file in euk_files:
            final_path = path + "/" + file
            df = pd.read_csv(final_path, sep = "\t")
            run = file.split("_")[0]
            output = generate_dicts(dict_taxid_sam, dict_taxid_lin, count_filt, count_read, count_success, count_na_rank, df, run, rank)
            count_filt, count_read, count_success , count_na_rank = output[0], output[1], output[2], output[3]
        
        #Show total files and messages for the directory-study
        tot_files = count_filt + count_read + count_success
        print("\n" + path)
        print("\nTotal de archivos procesados por EukDetect:", tot_files) 
        print("Archivos que no han conseguido superar los requisitos del filtro"
              " taxonomico:", count_filt)
        print("Archivos que presentan un recuento de lecturas vacio, probablemente porque no haya "
            "lecturas alineadas en la muestra:", count_read)
        print("Archivos que presentan una tabla vacia despues de seleccionar el Rank de interes:", count_na_rank)
        print("Archivos que se han procesado con exito:", count_success)
        
        
    #Generate final tables
    pres_abs_tab = pres_abs_table(dict_taxid_sam, dict_taxid_lin)
    taxonomy_tab = taxonomy_table(rank, dict_taxid_lin)


    #Name of the resulting tables
    if prefix_name:
        name_tab1 = prefix_name + "_presence_absence_" + rank + "_table.tsv"
        name_tab2 = prefix_name + "_taxonomy_" + rank + "_table.tsv"
    else:
        name_tab1 = "presence_absence_" + rank + "_table.tsv"
        name_tab2 = "taxonomy_" + rank + "_table.tsv"

    #Save the paths of result tables in a specified directory or in the currently one
    if output_dir:
        path_1 = os.path.join(output_dir, name_tab1)
        path_2 = os.path.join(output_dir, name_tab2)
    else:
        current_path = os.getcwd()
        path_1 = current_path + "/" + name_tab1
        path_2 = current_path + "/" + name_tab2

    #Save the results
    pres_abs_tab.to_csv(path_1, sep='\t',header=True, index=False)
    taxonomy_tab.to_csv(path_2, sep='\t',header=True, index_label="TaxID")

    tf = time()
    tot_time = tf - t0
    print("\nTiempo de ejecucion:", tot_time)
    print("\n")
    

if __name__ == '__main__':
    main()