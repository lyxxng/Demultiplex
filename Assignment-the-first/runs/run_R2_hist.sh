#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=R2
#SBATCH --output=R2-%j.out
#SBATCH --error=R2-%j.err
#SBATCH --mail-user=lxy@uoregon.edu
#SBATCH --mail-type=ALL

R2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
OUT="/projects/bgmp/lxy/bioinfo/Bi622/Demultiplex/Assignment-the-first/hists/R2"

/usr/bin/time -v python ../qual_dist.py -i $R2 -o $OUT -l 8 -n "Index 1"