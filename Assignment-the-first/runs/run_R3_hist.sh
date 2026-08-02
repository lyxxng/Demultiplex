#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=R3
#SBATCH --output=R3-%j.out
#SBATCH --error=R3-%j.err
#SBATCH --mail-user=lxy@uoregon.edu
#SBATCH --mail-type=ALL

R3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
OUT="/projects/bgmp/lxy/bioinfo/Bi622/Demultiplex/Assignment-the-first/hists/R3"

/usr/bin/time -v python ../qual_dist.py -i $R3 -o $OUT -l 8 -n "Index 2"