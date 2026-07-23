# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | phred+33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. 3,976,613 indexes in R2 and 3,328,051 indexes in R3 have undetermined base calls. This is 7,304,664 indexes in total.
    ```bash
    Commands used:
     $ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l
    3976613
    $ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l
    3328051
    ```
    
## Part 2
1. Define the problem

In multiplexing, multiple samples are pooled together so they can be sequenced in parallel on a single lane. This is done because it is efficient highly and cost-effective. Indexing can then be used to sort out reads after sequencing. After sequencing, the samples need to be linked back to their specific sample, or demultiplexed. This program takes the four reads (two sequences and two indexes), and sorts the sequences according to their indexes. It also handles index hopping, which is when a read has an incorrect index.

2. Describe output

* one R1 FASTQ file and one R2 FASTQ file per each matching index pair
* two FASTQ files for index hopping
* two FASTQ files for unknown (has at least one 'N') or low quality (does not meet quality score cutoff) reads
* tab-separated file with numbers of reads in each category, i.e.:
```
Total Matched   200
Total Hopped    100
Unknown 100
GTAGCGTA-GTAGCGTA  5
CTCTGGAT-ATCATGCG   0
...
```

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

**Answers for 4 and 5 located in [pseudocode.md](pseudocode.md)**