#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=demultiplex
#SBATCH --output=demultiplex-%j.out
#SBATCH --error=demultiplex-%j.err
#SBATCH --mail-user=lxy@uoregon.edu
#SBATCH --mail-type=ALL

DIR="/projects/bgmp/shared/2017_sequencing"
R1="$DIR/1294_S1_L008_R1_001.fastq.gz"
R2="$DIR/1294_S1_L008_R2_001.fastq.gz"
R3="$DIR/1294_S1_L008_R3_001.fastq.gz"
R4="$DIR/1294_S1_L008_R4_001.fastq.gz"
INDEX="$DIR/indexes.txt"
SCRATCH="/scratch/bgmp/lxy/demux"

/usr/bin/time -v python ./demultiplex.py -1 $R1 -2 $R2 -3 $R3 -4 $R4 -i $INDEX -o $SCRATCH -q 0