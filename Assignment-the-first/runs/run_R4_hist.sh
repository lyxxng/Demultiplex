#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=R4
#SBATCH --output=R4-%j.out
#SBATCH --error=R4-%j.err
#SBATCH --mail-user=lxy@uoregon.edu
#SBATCH --mail-type=ALL

R4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
OUT="/projects/bgmp/lxy/bioinfo/Bi622/Demultiplex/Assignment-the-first/hists/R4"

/usr/bin/time -v python ../qual_dist.py -i $R4 -o $OUT -l 101 -n "Read 2"