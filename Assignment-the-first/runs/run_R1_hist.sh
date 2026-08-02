#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=R1
#SBATCH --output=R1-%j.out
#SBATCH --error=R1-%j.err
#SBATCH --mail-user=lxy@uoregon.edu
#SBATCH --mail-type=ALL

R1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
OUT="/projects/bgmp/lxy/bioinfo/Bi622/Demultiplex/Assignment-the-first/hists/R1"

/usr/bin/time -v python ../qual_dist.py -i $R1 -o $OUT -l 101 -n "Read 1"