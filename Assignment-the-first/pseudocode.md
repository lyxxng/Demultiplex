pseudocode:
```
index set and dict function
    open and loop over indexes text file to read
        save indexes to a set
    create empty dict
    loop through itertools product of index set
        add indexes as key to dict as 'index1-index2' with value 0

unknown, matching, hopping num all initialized to 0

define dictionary for reverse compliment bp conversion (e.g., key A maps to value T, key T maps to value A)

reverse compliment function
    define an empty string
    loop over index in reverse
        append compliment bp to empty string
    return string

qual_score function from bioinfo module

demultiplex function
    open and loop over all four input fastq files for reading, output for writing
        if empty line
            break from while loop
        if on line with header
            save both headers in array
        if on line with seq
            convert index2 from reverse compliment
            save index1 and converted index2 in variables
            append index1-index2 to headers
        if on line with qual scores
            if 'N' in index1 or index2
                add to unknown files and increment unk num
            if index1 or index2 not in index set
                add to unknown files and increment unk num
            get average phred score
            if avg phred score below cut off
                add to unknown files and increment unk num
            if index1 is equal to index2
                add to correct matching files and increment match num
            else
                add to index hopping files and increment index hopping num

write results function
    open text file to write
        write matched, hopped, unk, and index dict numbers with tab separated
```
```python
def get_indexes(index_file: str) -> (set, dict):
    '''Takes a tab separate file with indexes and
    returns a set with each index and a dictionary with all pair combinations'''
    return index_set, index_dict
# Input: Text file with indexes "ATCG, GCAT"
# Expected output: ('ATCG', 'GCAT'), {'ATCG-GCAT': 0, 'ATCG-ATCG': 0, 'GCAT-GCAT': 0, 'GCAT-ATCG': 0}

def reverse_compliment(rev_seq: str) -> str:
    '''Takes a reverse compliment DNA string and
    returns the corresponding DNA sequence'''
    return seq
# Input: "ATGCTTGA"
# Expected output: "TCAAGCAT"

def demultiplex(qual_cutoff: int, R1: str, R2: str, R3: str, R4: str) -> (dict, int, int, int):
    '''Takes user defined quality score cutoff and input FASTA files and
    returns number of unknown, matched, and hopped'''
    return index_dict, matched, hopped, unk
# Input: qual_cutoff: 20, files in TEST-input_FASTQ
# Output: index_dict: {'GTAGCGTA-GTAGCGTA': 1, ...}, matched: 1, hopped: 1, unk: 1

def write_results(output_file: str, index_dict: dict, matched: int, hopped: int, unk: int):
    '''Takes output file name and numbers from demultiplex function and
    writes the information to a tab separate file'''
    return
# Input: index_dict: {'GTAGCGTA-GTAGCGTA': 5, 'CTCTGGAT-ATCATGCG': 2, ...}
#        matched: 5, hopped: 2, unk: 0
# Output: written to output_file tab separated text file --
# 
# Total Matched 5
# Total Hopped  2
# Unknown   0
# GTAGCGTA-GTAGCGTA  5
# CTCTGGAT-ATCATGCG   0
# Zero's for rest of pairs...
```