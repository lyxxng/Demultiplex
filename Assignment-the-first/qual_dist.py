#!/usr/bin/env python

import bioinfo
import gzip
import numpy as np

READ1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
INDEX1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
INDEX2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
READ2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

r1_arr, r2_arr = np.zeros(101), np.zeros(101)
i1_arr, i2_arr = np.zeros(8), np.zeros(8)

line = 0
with gzip.open(READ1, "rt") as r1, gzip.open(INDEX1) as i1, \
    gzip.open(INDEX2) as i2, gzip.open(READ2) as r2:
    while True:
        r1_line = r1.readline().strip("\n")
        if r1_line == "":
            break
        # print(r1_line)
        if line % 4 == 3:
            # print(r1_line)
            for nt, qual in enumerate(r1_line):
                # print(str(qual))
                r1_arr[nt] += bioinfo.convert_phred(qual)
        line += 1

print(r1_arr)