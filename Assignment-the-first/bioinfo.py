#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set("ATGCNatgcn")
RNA_bases = set("AUGCNaugcn")

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    '''Computes the average phred quality score given a string of phred scores'''
    total = 0
    length = len(phred_score)

    for i in range(length):
        total += convert_phred(phred_score[i])

    return total / length

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq) <= (RNA_bases if RNAflag else DNA_bases)

def gc_content(seq: str, RNAflag: bool=False) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq, RNAflag), "String contains invalid characters"

    return (seq.count("G") + seq.count("C")) / len(seq)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    num_values = len(lst)
    middle = num_values // 2

    if num_values % 2 == 0:
        return (lst[middle] + lst[middle - 1]) / 2
    else:
        return lst[middle]

def oneline_fasta(in_file, parsed_file):
    '''Convert FASTA file so it has 2 lines per record (combine seq data into one line)'''
    with open(in_file, "r") as in_f, open(parsed_file, "w") as out_f:
        first_line = True

        for line in in_f:
            if line[0] == ">":
                # Prepend newline if not the very first header
                if first_line:
                    first_line = False
                    out_f.write(line)
                else:
                    out_f.write(f"\n{line}")
            else:
                out_f.write(line.strip("\n"))
        out_f.write("\n")

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert qual_score("ABC") == 33.0
    assert qual_score("#*+5:") == 13.2
    assert qual_score("HIJK") == 40.5
    assert qual_score("GGG??") == 34.8
    print("Your qual_score function is working!")

    assert validate_base_seq("ATAGCGT")
    assert validate_base_seq("AGUCUCUCUAG", True)
    assert validate_base_seq("augucgcu", True)
    assert validate_base_seq("ksdjfw") == False
    print("Your validate_base_seq function is working!")

    assert gc_content("GCGC") == 1
    assert gc_content("AUAUAUAU", True) == 0
    assert gc_content("GACAUUCAUAAU", True) == 0.25
    assert gc_content("GCGCGCAT") == 0.75
    print("Your gc_content function is working!")

    assert calc_median([9, 10, 12, 13, 13, 13, 15, 15, 16, 16, 18, 22, 23, 24, 24, 25]) == 15.5
    assert calc_median([4, 40, 80, 87, 96]) == 80
    assert calc_median([0, 18, 26]) == 18
    assert calc_median([13, 95]) == 54
    assert calc_median([4, 8, 22, 46, 53, 71]) == 34
    print("Your calc_median function is working!")