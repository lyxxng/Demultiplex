#### 07/21
Initial data exploration. Index1 match given indexes, Index2 is the reverse compliment:
```bash
$ zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep -v "^N" | head
TACCGGAT
CTCTGGAT
CTCTGGAT
TACCGGAT
CTCTGGAT
TACCGGAT
TACCGGAT
AGAGTCCA
TCTTCGAC
GTAGCGTA

$ zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep -v "^N" | head
GCTATCCT
ATCCAGAG
TGAGCTAG
CATGGCCG
ATCCAGAG
TGAGCTAG
ATCGATCG
ATCCAGAG
TGGACTCT
ACTCTCGA
```
Records appear to come in the same order for read files and index files. Also, '#' encodes for 'N', implying Phred+33 encoding:
```bash
$ zcat 1294_S1_L008_R1_001.fastq.gz | head
@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ
@K00337:83:HJKJNBBXX:8:1101:1286:1191 1:N:0:1
CNACCTGTCCCCAGCTCACAGGACAGCACACCAAAGGCGGCAACCCACACCCAGTTTTACAGCCACACAGTGCCTTGTTTTACTTGAGGACCCCCCACTCC
+
A#AAFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJAJJJJJJJJJJJJJJFJJJJJFFFFJJJJJJJJJJJJJJJJJJ77F
@K00337:83:HJKJNBBXX:8:1101:1347:1191 1:N:0:1
GNGGTCTTCTACCTTTCTCTTCTTTTTTGGAGGAGTAGAATGTTGAGAGTCAGCAGTAGCCTCATCATCACTAGATGGCATTTCTTCTGAGCAAAACAGGT

$ zcat 1294_S1_L008_R2_001.fastq.gz | head
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ
@K00337:83:HJKJNBBXX:8:1101:1286:1191 2:N:0:1
NACAGCGA
+
#AAAFJJJ
@K00337:83:HJKJNBBXX:8:1101:1347:1191 2:N:0:1
NTCCTAAG
```
Read length is 101 and 8 for reads and indexes, respectively:
```bash
$ zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc -c
102

$ zcat 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc -c
9
```
363,246,735 records:
```bash
$ zcat 1294_S1_L008_R1_001.fastq.gz | wc -l
1452986940

1452986940 / 4 = 363246735
```

#### 07/22
Finding how many indexes have undetermined base calls:
```bash
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l
3976613

$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l
3328051
```

Manually created test files located in `TEST-input_FASTQ/` and `TEST-output_FASTQ`

#### 08/01
Wrote a script for making the quality score distributions located at `Assignment-the-fist/qual_dist.py`. I tested it with a smaller file first by piping the first 100 records:
```bash
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -400 > Assignment-the-first/test.fq
```
Once I determined that my script worked correctly for the small test file, I ran the script on the four larger files. I used argparse so I could input the name of the input and output files and some other data used in the code, and ran the script four times as batch jobs. Here is the `usr/bin/time` data from those runs:
| Job | Elapsed time | CPU | MRSS |
| --- | --- | --- | --- |
| run_R1_hist.sh | 1:44:41 | 99% | 70.4 MB |
| run_R2_hist.sh | 13:11.16 | 99% | 72.3 MB |
| run_R3_hist.sh | 13:13.73 | 99% | 80.0 MB |
| run_R4_hist.sh | 1:45:31 | 99% | 70.2 MB |

#### 08/04
Ran my demultiplexing script. Job number was 45990590.
| Job | Elapsed time | CPU | MRSS |
| --- | --- | --- | --- |
| run_demultiplex.sh | 1:08:32 | 75% | 258.4 MB |

Used a quality cutoff of 30 and had a lot of unknowns. Going to rewrite to use no quality cutoff and also use lists instead of numpy arrays, since I think that is slowing it down.

#### 08/05
Rewrote my code to store the records in lists instead of numpy arrays. Job number was 46007752. It was 20 minutes faster than the previous version.
| Job | Elapsed time | CPU | MRSS |
| --- | --- | --- | --- |
| run_demultiplex.sh | 48:17.46 | 78% | 246.0 MB |