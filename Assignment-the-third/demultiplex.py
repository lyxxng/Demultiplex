#!/usr/bin/env python

import itertools

index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Add indexes to a set and all index pair combinations to a dictionary
def get_indexes(index_file):
    index_set = set()

    with open(index_file, "r") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            index_set.add(line.strip("\n").split("\t")[4])

    index_dict = {}
    for pair in itertools.product(index_set, repeat=2):
        index_dict[pair[0] + "-" + pair[1]] = 0
    # print(index_dict)

    return index_set, index_dict

# Get reverse compliment of a given index/barcode
def reverse_compliment(index):
    rev_comp = ""
    for bp in index[::-1]:
        rev_comp += comp_dict[bp]

    return rev_comp

# print(rev_comp)