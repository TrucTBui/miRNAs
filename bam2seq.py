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
            chromosome = columns[2]
            start_position = int(columns[3])  # Start position (1-based)
            sequence = columns[9]          # Sequence of the read

            reads.append([chromosome, start_position, sequence])

    return reads


def trim_reads_and_find_variants(all_reads, chromosome, start_range, end_range):
    variant_dict = {}

    for read in all_reads:
        read_chromosome = read[0]
        read_start = read[1]  # Start pos of read
        read_end = read_start + len(read[2]) - 1

        # check for suitable reads
        if chromosome != read_chromosome:
            continue
        if read_end < start_range or read_start > end_range:
            continue

        # determine the actual start and end positions of the trimmed read
        actual_start = max(read_start, start_range)
        actual_end = min(read_end, end_range)

        # trim the sequence
        trimmed_sequence = read[2][actual_start - read_start:actual_end - read_start + 1]

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
        frequency = Counter(bases)  # Count how frequent each variant is
        unique_bases = set(bases)
        if len(unique_bases) > 1:  # More than one unique nucleotide means a variant
            sorted_bases = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
            ordered_bases = [base for base, count in sorted_bases]
            seq += "(" + "/".join(ordered_bases) + ")"
            variants[pos] = dict(sorted_bases)
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
    genome_start_time = time.perf_counter()  # runtime measurement

    # Extract all reads for the specified locations
    combined_locations = ' '.join(locations)
    temp_sam_path = f"/tmp/{genome_basename}.sam"
    command1 = f"(cd /mnt/proj/software && samtools view {genome} {combined_locations}) > {temp_sam_path}"
    subprocess.run(command1, shell=True, check=True)

    all_reads = extract_reads_from_sam(temp_sam_path)

    for location in locations:
        # Parse chromosome and range
        chrom, pos_range = location.split(":")
        start_pos, end_pos = map(int, pos_range.split("-"))

        # Process reads and find variants for the current location
        seq, variants = trim_reads_and_find_variants(all_reads, chrom, start_pos, end_pos)

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
    os.remove(f"/tmp/{genome_basename}.sam")
    genome_end_time = time.perf_counter()
    runtime = genome_end_time - genome_start_time
    print(f"Runtime for {genome_basename}: {runtime:.5f} seconds")

# Write to output file if specified
if output:
    with open(output, "w") as o:
        o.write("\n".join(results))
