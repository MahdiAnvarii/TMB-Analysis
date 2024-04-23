#!/usr/bin/env python3

import sys

def count_sequencing_coverage(bed_file):
    """Function to count the total number of sequencing coverage from a bed file."""
    coverage = 0
    with open(bed_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            coverage += int(fields[2]) - int(fields[1])
    return coverage

def calculate_TMB(total_mutations):
    """Function to calculate Tumor Mutational Burden."""
    coverage = count_sequencing_coverage('/mnt/d/NGS/Samples/P2/Data/SRR28000175_sorted_dedup_bqsr_merged.bed')
    # print(coverage)
    tmb = total_mutations / (coverage / 1000000)
    return tmb

def main():
    if len(sys.argv) != 2:
        print("Usage: {} <total_mutations>".format(sys.argv[0]))
        sys.exit(1)

    total_mutations = int(sys.argv[1])
    tmb = calculate_TMB(total_mutations)
    print("Tumor Mutational Burden (TMB): {:.2f} mutations per megabase".format(tmb))

if __name__ == "__main__":
    main()
