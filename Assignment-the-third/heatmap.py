#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns

index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

indexes = []
with open(index_file, "r") as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        indexes.append(line.strip("\n").split("\t")[4])

counts = np.empty([24, 24])

with open("results.tsv", "r") as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.strip("\n").split("\t")

        index1, index2 = line[0].split("-")

        counts[indexes.index(index1), indexes.index(index2)] = int(line[1]) / 363246735 * 100

fig, ax = plt.subplots()

sns.heatmap(counts, xticklabels=indexes, yticklabels=indexes,
            norm=mcolors.LogNorm(), cmap="coolwarm", linewidths=0.5)

plt.ylabel("Index 1")
plt.xlabel("Index 2")
ax.set_title("Heatmap of Index Pairs")
fig.tight_layout()
plt.savefig("heatmap.png")