#!/usr/bin/env python

import argparse
import bioinfo
import gzip
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="name of file for input")
    parser.add_argument("-o", help="name of file for output")
    parser.add_argument("-l", help="length of sequence")
    parser.add_argument("-n", help="name of read (e.g., Read 1 or Index 1)")

    return parser.parse_args()

args = get_args()

qual_arr = np.zeros(int(args.l))

line = 0
with gzip.open(args.i, "rt") as f:
    for i, line in enumerate(f):
        # Grab only the lines with quality score info
        if i % 4 == 3:
            line = line.strip("\n")
            for nt, qual in enumerate(line):
                qual_arr[nt] += bioinfo.convert_phred(qual)
    # Each index contains the average qual score
    qual_arr = qual_arr / (i / 4)

plt.scatter(range(int(args.l)), qual_arr)
plt.title(f"Quality Score Distribution Per Nucleotide for {args.n}")
plt.xlabel("Nucleotide Index")
plt.ylabel("Average Phred Quality Score")
plt.savefig(f"{args.o}_hist.png")