#!/usr/bin/env python

import bioinfo
import gzip
import itertools
import numpy as np

index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
DIR="/projects/bgmp/shared/2017_sequencing/"
DIR="/projects/bgmp/lxy/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/"
SCRATCH = "/scratch/bgmp/lxy/demux/"

# Complement of each base
comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

# TODO: argparse()
# [ ] all four fq input files
# [ ] index file
# [ ] output directory (i.e., scratch)
# [ ] quality score cutoff

# TODO: write results function

# Open all necessary files
def open_files(index_set):
    file_dict = {}

    # Open all matching files to write
    for i in index_set:
        file_dict[i] = [open(f"{SCRATCH}{i}.R1.fastq", "w")]
        file_dict[i].append(open(f"{SCRATCH}{i}.R2.fastq", "w"))

    # Open all index hopped files to write
    file_dict["hopped"] = [open(f"{SCRATCH}hopped.R1.fastq", "w")]
    file_dict["hopped"].append(open(f"{SCRATCH}hopped.R2.fastq", "w"))

    # Open all unknown index files to write
    file_dict["unk"] = [open(f"{SCRATCH}unknown.R1.fastq", "w")]
    file_dict["unk"].append(open(f"{SCRATCH}unknown.R2.fastq", "w"))

    ###################### TEST FILES, REMOVE ######################
    file_dict["help"] = [open(f"help.R1.fastq", "w")]
    file_dict["help"].append(open(f"help.R2.fastq", "w"))

    file_dict["match"] = [open(f"match.R1.fastq", "w")]
    file_dict["match"].append(open(f"match.R2.fastq", "w"))

    file_dict["hop"] = [open(f"hop.R1.fastq", "w")]
    file_dict["hop"].append(open(f"hop.R2.fastq", "w"))

    # read1 = gzip.open(f"{DIR}1294_S1_L008_R1_001.fastq.gz", "rt")
    # index1 = gzip.open(f"{DIR}1294_S1_L008_R2_001.fastq.gz", "rt")
    # index2 = gzip.open(f"{DIR}1294_S1_L008_R3_001.fastq.gz", "rt")
    # read2 = gzip.open(f"{DIR}1294_S1_L008_R4_001.fastq.gz", "rt")

    # Open all input files to read
    read1 = open(f"{DIR}read1.fastq", "r")
    index1 = open(f"{DIR}index1.fastq", "r")
    index2 = open(f"{DIR}index2.fastq", "r")
    read2 = open(f"{DIR}read2.fastq", "r")

    # Return dictionary of write files and all input files as a list
    return file_dict, [read1, index1, index2, read2]

# Close all files at the end of program
def close_files(file_dict, in_files):
    # Output files
    for key in file_dict:
        file_dict[key][0].close()
        file_dict[key][1].close()

    # Input files
    for file in in_files:
        file.close()

    return

# Add indexes to a set and all index pair combinations to a dictionary
def get_indexes(index_file):
    index_set = set()

    with open(index_file, "r") as f:
        for i, line in enumerate(f):
            # Skip header line
            if i == 0:
                continue
            # Get just the index string
            index_set.add(line.strip("\n").split("\t")[4])

    # Save all possible combinations of index pairs
    index_dict = {}
    for pair in itertools.product(index_set, repeat=2):
        key = f"{pair[0]}-{pair[1]}"
        index_dict[key] = 0

    return index_set, index_dict

# Get reverse compliment of a given index/barcode
def reverse_compliment(index):
    rev_comp = ""
    for bp in index[::-1]:
        rev_comp += comp_dict[bp]

    return rev_comp

def demultiplex(index_set, index_dict, file_dict, in_files, qual_cutoff):
    matched, hopped, unknown = 0, 0, 0
    records = np.empty(shape=(4, 4), dtype=np.dtypes.StringDType)
    header = ""

    while True:
        # Save one record at a time for each file
        for i, file in enumerate(in_files):
            for j in range(4):
                records[i, j] = file.readline().strip("\n")

                # Exit loop if we reach end of file
                if i == 0 and j == 0 and records[i, j] == "":
                    return matched, hopped, unknown
                
        records[2, 1] = reverse_compliment(records[2, 1])

        # Add sequence of "index1-index2" to header of both reads
        header = f"{records[1, 1]}-{records[2, 1]}"
        records[0, 0] += f" {header}"
        records[3, 0] += f" {header}"

        # Write records to unknown if index does not exist, this should also handle if there is an 'N'
        if records[1, 1] not in index_set or records[2, 1] not in index_set:
            unknown += 1
            for i in range(4):
                file_dict["help"][0].write(f"{str(records[0, i])}\n")
                file_dict["help"][1].write(f"{str(records[3, i])}\n")
        # Write records to unknown if either average index quality score is below given cutoff
        elif bioinfo.qual_score(records[1, 3]) < qual_cutoff or bioinfo.qual_score(records[2, 3]) < qual_cutoff:
            unknown += 1
            for i in range(4):
                file_dict["help"][0].write(f"{str(records[0, i])}\n")
                file_dict["help"][1].write(f"{str(records[3, i])}\n")
        # Write records to matched if the indexes are the same
        elif records[1, 1] == records[2, 1]:
            matched += 1
            index_dict[header] += 1
            for i in range(4):
                file_dict["match"][0].write(f"{str(records[0, i])}\n")
                file_dict["match"][1].write(f"{str(records[3, i])}\n")
        # Write records to hopped if the indexes are not the same
        else:
            hopped += 1
            index_dict[header] += 1
            for i in range(4):
                file_dict["hop"][0].write(f"{str(records[0, i])}\n")
                file_dict["hop"][1].write(f"{str(records[3, i])}\n")

def main():
    # TODO: Call argparse()

    index_set, index_dict = get_indexes(index_file)

    file_dict, in_files = open_files(index_set)
    demultiplex(index_set, index_dict, file_dict, in_files, 30)
    close_files(file_dict, in_files)

    # TODO: Write relevant informaton to file

if __name__ == "__main__":
    main()