"""
tmp version 31.10.2024
Input: Genome file(s) (.bam), searched region(s)
Output: Reads within the given region with their variants
"""
import subprocess
import argparse
import os
from collections import Counter
import time


def extract_reads_from_sam(sam_file):
    reads = []

    with open(sam_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            start_position = int(columns[3])  # Start position (1-based)
            sequence = columns[9]          # Sequence of the read

            reads.append([start_position, sequence])

    return reads


def trim_reads_and_find_variants(all_reads, start_range, end_range):
    variant_dict = {}

    for read in all_reads:
        read_start = read[0]  # Start pos of read
        read_end = read_start + len(read[1]) - 1

        # double check
        if read_end < start_range or read_start > end_range:
            continue

        # determine the actual start and end positions of the trimmed read
        actual_start = max(read_start, start_range)
        actual_end = min(read_end, end_range)

        # trim the sequence
        trimmed_sequence = read[1][actual_start - read_start:actual_end - read_start + 1]

        # find variants
        for i, base in enumerate(trimmed_sequence):
            pos = actual_start + i
            if pos not in variant_dict:
                variant_dict[pos] = []
            variant_dict[pos].append(base)

    seq = ""
    variants = {}
    for pos in sorted(variant_dict.keys()):
        bases = variant_dict[pos]
        unique_bases = set(bases)
        if len(unique_bases) > 1:  # More than one unique nucleotide means a variant
            seq += "(" + "/".join(unique_bases) + ")"

            frequency = Counter(bases)  # count how frequent each variant is
            variants[pos] = dict(frequency)
        else:
            seq += bases[0]

    return seq, variants


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", type=str, nargs="+", required=True, help="genome file in bam format")
parser.add_argument("-l", "--location", type=str, nargs="+", required=True, help="genome location to extract reads, e.g. 9:1000-1005 for positions 1000-1005 of chromosome 9")
parser.add_argument("-o", "--output", type=str, required=False, help="output file for the result")

args = parser.parse_args()
genomes = args.genome
locations = args.location
output = args.output

results = [f"Genome\tLocation\tSequence\tVariants"]

for genome in genomes:
    genome_basename = os.path.splitext(os.path.basename(genome))[0]
    genome_start_time = time.time() #runtime measurement

    for location in locations:
        # Parse chromosome and range
        chrom, pos_range = location.split(":")
        start_pos, end_pos = map(int, pos_range.split("-"))

        # Temporary SAM file path for each genome-location combination
        temp_sam_path = f"/tmp/{genome_basename}_{chrom}_{start_pos}_{end_pos}.sam"

        # Extract reads from the BAM file using samtools
        command1 = f"(cd /mnt/proj/software && samtools view {genome} '{location}') > {temp_sam_path}"
        subprocess.run(command1, shell=True, check=True)

        # Process reads and find variants
        reads = extract_reads_from_sam(temp_sam_path)
        seq, variants = trim_reads_and_find_variants(reads, start_pos, end_pos)

        genome_end_time = time.time()
        runtime = genome_end_time - genome_start_time
        print(f"Runtime for {genome_basename} at {location}: {runtime:.2f} seconds")

        # Store results in structured format
        if output:
            results.append(f"{genome_basename}\t{location}\t{seq}\t" +
                           "; ".join(f"{pos}:{counts}" for pos, counts in variants.items()))
        else:
            print(f"Genome file: {genome}\nLocation: {location}")
            print(f"Sequence: {seq}")
            print("Variants and frequencies:")
            for pos, counts in variants.items():
                print(f"Position {pos}: {counts}")

        # Clean up
        os.remove(temp_sam_path)

# Write to output file if specified
if output:
    with open(output, "w") as o:
        o.write("\n".join(results))
