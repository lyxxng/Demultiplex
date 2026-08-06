#!/usr/bin/env python

import argparse
import bioinfo
import gzip
import itertools

# Complement of each base
comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def get_args():
    parser = argparse.ArgumentParser(description="A script to demultiplex sequencing data")

    parser.add_argument("-1", "--run1", help="input file for run 1 (read 1)", required=True)
    parser.add_argument("-2", "--run2", help="input file for run 2 (index 1)", required=True)
    parser.add_argument("-3", "--run3", help="input file for run 3 (index 2)", required=True)
    parser.add_argument("-4", "--run4", help="input file for run 4 (read 2)", required=True)
    parser.add_argument("-i", "--indexes", help="tsv file with indexes", required=True)
    parser.add_argument("-o", "--output_dir", help="path to output directory", default=".", required=False)
    parser.add_argument("-q", "--qual_cutoff", help="quality score cutoff", default=0, required=False)

    return parser.parse_args()

# Add indexes to a set and all index pair combinations to a dictionary
def get_indexes(index_file):
    index_set = set()

    with open(index_file, "r") as f:
        for i, line in enumerate(f):
            # Skip header line
            if i == 0:
                continue
            # Get just the index string
            index = line.strip("\n").split("\t")[4]
            index_set.add(index)

    # Save all possible combinations of index pairs
    index_dict = {}
    for pair in itertools.product(index_set, repeat=2):
        key = f"{pair[0]}-{pair[1]}"
        index_dict[key] = 0

    return index_set, index_dict

# Open all necessary files
def open_files(index_set, output_dir, r1, r2, r3, r4):
    file_dict = {}

    # Open all matching files to write
    for i in index_set:
        file_dict[i] = [open(f"{output_dir}/{i}.R1.fastq", "w")]
        file_dict[i].append(open(f"{output_dir}/{i}.R2.fastq", "w"))

    # Open all index hopped files to write
    file_dict["hopped"] = [open(f"{output_dir}/hopped.R1.fastq", "w")]
    file_dict["hopped"].append(open(f"{output_dir}/hopped.R2.fastq", "w"))

    # Open all unknown index files to write
    file_dict["unk"] = [open(f"{output_dir}/unknown.R1.fastq", "w")]
    file_dict["unk"].append(open(f"{output_dir}/unknown.R2.fastq", "w"))

    # Open all input files to read
    read1 = gzip.open(r1, "rt")
    index1 = gzip.open(r2, "rt")
    index2 = gzip.open(r3, "rt")
    read2 = gzip.open(r4, "rt")

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

# Get reverse compliment of a given index/barcode
def reverse_compliment(index):
    rev_comp = ""
    for bp in index[::-1]:
        rev_comp += comp_dict[bp]

    return rev_comp

def demultiplex(index_set, index_dict, file_dict, in_files, qual_cutoff):
    matched, hopped, unknown = 0, 0, 0
    header = ""

    while True:
        records = [[] for _ in range(4)]
        # Save one record at a time for each file
        for i, file in enumerate(in_files):
            for _ in range(4):
                line = file.readline().strip()
                records[i].append(line)

                # Exit loop if we reach end of file
                if line == "":
                    return matched, hopped, unknown
                
        records[2][1] = reverse_compliment(records[2][1])

        # Add sequence of "index1-index2" to header of both reads
        header = f"{records[1][1]}-{records[2][1]}"
        records[0][0] += f" {header}"
        records[3][0] += f" {header}"

        # Write records to unknown if index does not exist, this should also handle if there is an 'N'
        if records[1][1] not in index_set or records[2][1] not in index_set:
            unknown += 1

            file_dict["unk"][0].writelines(line + '\n' for line in records[0])
            file_dict["unk"][1].writelines(line + '\n' for line in records[3])
        # Write records to unknown if either average index quality score is below given cutoff
        elif qual_cutoff:
            if bioinfo.qual_score(records[1][3]) < qual_cutoff or bioinfo.qual_score(records[2][3]) < qual_cutoff:
                unknown += 1
                file_dict["unk"][0].writelines(line + '\n' for line in records[0])
                file_dict["unk"][1].writelines(line + '\n' for line in records[3])
        # Write records to matched if the indexes are the same
        elif records[1][1] == records[2][1]:
            matched += 1
            index_dict[header] += 1
            file_dict[records[1][1]][0].writelines(line + '\n' for line in records[0])
            file_dict[records[1][1]][1].writelines(line + '\n' for line in records[3])
        # Write records to hopped if the indexes are not the same
        else:
            hopped += 1
            index_dict[header] += 1
            file_dict["hopped"][0].writelines(line + '\n' for line in records[0])
            file_dict["hopped"][1].writelines(line + '\n' for line in records[3])

def write_results(index_dict, matched, hopped, unknown):
    num_records = matched + hopped + unknown

    with open("results.md", "w") as md, open ("results.tsv", "w") as tsv:
        md.write("### Main stats\n| | Total | Percentage |\n| --- | --- | --- |\n")
        md.write(f"| **Matched** | {matched} | {round(matched / num_records * 100, 2)} |\n")
        md.write(f"| **Hopped** | {hopped} | {round(hopped / num_records * 100, 2)} |\n")
        md.write(f"| **Unknown** | {unknown} | {round(unknown / num_records * 100, 2)} |\n")

        md.write("---\n### Stats per barcode pair\nNote: Matching barcodes are in bold.\n")
        md.write("| | Total | Percentage |\n| --- | --- | --- |\n")
        
        tsv.write("Barcode\tCount\n")
        for key, value in index_dict.items():
            tsv.write(f"{key}\t{value}\n")

            indexes = key.split("-")
            if indexes[0] == indexes[1]:
                md.write(f"| **{key}** | **{value}** | **{round(value / num_records * 100, 2)}** | \n")
            else:
                md.write(f"| {key} | {value} | {round(value / num_records * 100, 2)} | \n")

        return

def main():
    args = get_args()

    r1, r2, r3, r4 = args.run1, args.run2, args.run3, args.run4
    index_file, output_dir = args.indexes, args.output_dir
    qual_cutoff = int(args.qual_cutoff)

    index_set, index_dict = get_indexes(index_file)

    file_dict, in_files = open_files(index_set, output_dir, r1, r2, r3, r4)
    matched, hopped, unknown = demultiplex(index_set, index_dict, file_dict, in_files, qual_cutoff)
    close_files(file_dict, in_files)

    write_results(index_dict, matched, hopped, unknown)

if __name__ == "__main__":
    main()