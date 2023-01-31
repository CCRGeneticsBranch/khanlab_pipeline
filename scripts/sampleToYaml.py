#!/usr/bin/python
import re
import requests
import os,subprocess
import pandas as pd
import yaml
import sys
import argparse
import re

def main(args):

    #User input:
    #1. sample_id
    
    sample_id = args.sample.strip()
    out_file = args.out.strip()
    masterfile_dir="/data/khanlab/projects/DATA/Sequencing_Tracking_Master/"
    master_files = ["Sequencing_Tracking_Master_db.txt","ClinOmics_Sequencing_Master_File_db.txt","SequencingMasterFile_OutsidePatients_db.txt"]
    #sample information from Young's master file
    sample_file = "Sample_" + sample_id
    master_sample = {}
    default_genome="hg19"
    xeno_genome="mm10"
    #columns = ["Type of sequencing","Matched normal","Matched RNA-seq"]
    
    for master_file in master_files:        
        file = masterfile_dir + "/" + master_file
        master_df = pd.read_csv(file, sep='\t', encoding = "ISO-8859-1")
        master_df['Sample_ID'] = master_df['Library ID'] + '_' + master_df['FCID']
        master_df['SampleFiles'] = 'Sample_' + master_df['Library ID'] + '_' + master_df['FCID']
        sample_df = master_df.loc[(master_df['Sample_ID'] == sample_id) | (master_df['Library ID'] == sample_id)]
        sample_df = sample_df.reset_index(drop=True)
        if sample_df['Sample_ID'].count() or sample_df['Library ID'].count() > 0:            
            for column in sample_df:
                master_sample[column] = str(sample_df[column][0]).strip()
            if master_sample["Sample_ID"] == "nan":
                master_sample["Sample_ID"] = master_sample["Library ID"]
                master_sample["SampleFiles"] = 'Sample_' + master_sample["Library ID"]
            matched_normal = master_sample["Matched normal"]
            matched_df = master_df.loc[(master_df['Sample_ID'] == matched_normal) | (master_df['Library ID'] == matched_normal)]
            matched_df = matched_df.reset_index(drop=True)
            if matched_df['Sample_ID'].count() > 0:
                master_sample["Matched normal"] = str(matched_df["Sample_ID"][0])
            matched_rnaseq = master_sample["Matched RNA-seq lib"]
            matched_rnaseq_df = master_df.loc[(master_df['Sample_ID'] == matched_rnaseq) | (master_df['Library ID'] == matched_rnaseq)]
            matched_rnaseq_df = matched_rnaseq_df.reset_index(drop=True)
            if matched_rnaseq_df['Sample_ID'].count() > 0:
                master_sample["Matched RNA-seq lib"] = str(matched_rnaseq_df["Sample_ID"][0])
            if master_sample["Type"].find("xeno") >=0:
                master_sample["Xenograft"] = "yes"
                master_sample["XenograftGenome"] = xeno_genome
            if master_sample["Type of sequencing"] == "C-il":
                print("chipseq")
            elif master_sample["Type of sequencing"] == "H-il":
                print("hic")
            elif "T-il" in master_sample["Type of sequencing"]:
                print("rnaseq")
            else:
                print("Sequencing type " + master_sample["Type of sequencing"] + " is not supported")
                exit(1)
            if master_sample["SampleRef"] == "nan":
                master_sample["SampleRef"] = "hg19"            
            print(master_sample["SampleRef"])
            break
    if 'Sample_ID' not in master_sample.keys():
        print(sample_id + " not found in Khanlab master files\n")
        sys.exit(1)
    fo = open(out_file,"w")
    fo.write(yaml.dump({"samples":{master_sample["Sample_ID"] : master_sample}}))
    fo.close()
            
    
parser = argparse.ArgumentParser(description='Generate YAML sample sheet.')
parser.add_argument("--sample", "-s", metavar="SAMPLE_ID", required=True, help="Sample ID (Library_id)_(FCID)")
parser.add_argument("--out", "-o", metavar="OUTPUT", required=True, help="Output YAML file")
#parser.add_argument("--chip", "-c", metavar="ChIPseq MASTER", help="Chipseq Master file", default="/data/khanlab/projects/ChIP_seq/manage_samples/ChIP_seq_samples.xlsx")
#parser.add_argument("--hic", "-i", metavar="HiC MASTER", help="HiC Master file", default="/data/khanlab/projects/HiC/manage_samples/HiC_master.xlsx")
args = parser.parse_args()

main(args)