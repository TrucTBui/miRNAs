"""
version 22.11.2024
Input: A file containing paths to genomes, A file containing paths to searched region(s), Path to reference genome
Output: Reads within the given region with their variants in sequencing and in comparison with the reference genome
"""
import subprocess
import argparse
import os
from collections import Counter
import time
import pandas as pd


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

def reads_processing(all_reads, location):
    chromosome, pos_range = location.split(":")
    start_range, end_range = map(int, pos_range.split("-"))

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
            # pos = f"{chromosome}:{actual_start + i}"
            pos = actual_start + i
            if pos not in variant_dict:
                variant_dict[pos] = []
            variant_dict[pos].append(base)

    return variant_dict

def find_variants_per_pos(variant_dict: dict(), reference):
    seq = ""
    variants = {}  # Contains base of reference and alternative
    all_results = {}
    SNPs = {}
    possible_SNPs = {}
    counter = 0  # for iterating the seq from ref genome

    for pos in sorted(variant_dict.keys()):
        ref = reference[counter]
        all_results[pos] = [ref]

        bases = variant_dict[pos]
        frequency = Counter(bases)  # Count how frequent each variant is
        unique_bases = set(bases)
        all_results[pos].append("/".join(unique_bases))

        # Find variances in reads
        if len(unique_bases) > 1:  # More than one unique nucleotide means a variant
            sorted_bases = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
            ordered_bases = [base for base, count in sorted_bases]
            seq += "(" + "/".join(ordered_bases) + ")"
            variants[pos] = dict(sorted_bases)

            # Check if it can be a SNP or Haploid difference
            most_frequent_base = ordered_bases[0]  # Extract the most significant base
            ratio = (frequency[most_frequent_base] / sum(frequency.values()))

            # compare the most frequent base read at position with the ref
            if not ordered_bases[0] == ref or ratio < 0.7:
                all_bases = "/".join(ordered_bases)
                possible_SNPs[pos] = f"{ref}//{all_bases}"

        else:  # Only one base at that position
            base = bases[0]
            seq += base

            if not base == ref:  # SNP found
                SNPs[pos] = f"{ref}/{base}"

        counter += 1

    return all_results, seq, variants, SNPs, possible_SNPs

def extract_ref_genome(fa, locations):
    """
    :param fa: path to fasta file containing reference genome
    :param locations: a file containing loctions to extract
    :return: a dictionary, where key=location, value=sequence at the location
    """
    cmd = f"samtools faidx {fa} -r {locations}"
    result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
    output = result.stdout.decode("utf-8")  # Decode the `stdout` attribute

    seqs= {}
    seq = ""
    location = None

    for line in output.splitlines():
        if line.startswith(">"):  # Header line
            if location:
                seqs[location] = seq
            location = line[1:]  # Extract location name (without '>')
            seq = ""
        else:
            seq += line.strip()

    # Add the last sequence to the dictionary
    if location:
        seqs[location] = seq

    return seqs


def parse_genome(path):
    genomes = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            genomes.append(line.strip())
    return genomes

def calc_false_sequencing_rate(variants:dict, possible_SNPs:dict, ref):
    if len(variants) > 0:
        return round((len(variants) - len(possible_SNPs)) / len(ref), 5)
    else:
        return 0.0

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", type=str, required=True, help="file containing path to genome(s) in bam format")
parser.add_argument("-r", "--reference", type=str, required=True, help="reference genome in fasta format")
parser.add_argument("-l", "--location", type=str, required=True, help="genome location to extract reads")
parser.add_argument("-o", "--output", type=str, required=False, help="output folder for the result", default="./")

args = parser.parse_args()
genome_path = args.genome
ref_path = args.reference
output = args.output
location_path = args.location

# Parse the input path to a list of bam files (index 0) and vcf files(index1)
genomes= parse_genome(genome_path)

locations=[]

if os.path.isfile(location_path):  # Input is a file
    with open(location_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                locations.append(line.strip())
else:
    raise Exception("Location file not found or not correct")

ref_genome = extract_ref_genome(ref_path, location_path)

results = [f"Genome\tChromosome\tPosition\tReference\tAlternative"]

summary = [f"Genome\tChromosome\tStart\tEnd\tSequence\tReference\tSNP\tSNP_Freq\tRead_Variance\tVariance_Freq\tPossible_SNP\tApprox_sequencing_rate"]

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
        # Extract ref genome at the current range
        ref = ref_genome.get(location)

        # Process reads and find variants for the current location
        variant_dict = reads_processing(all_reads, location)

        all_results, seq, variants, SNPs, possible_SNPs = find_variants_per_pos(variant_dict, ref)

        chrom, pos_range = location.split(":")
        start_pos, end_pos = map(int, pos_range.split("-"))

        false_seq_rate = calc_false_sequencing_rate(variants, possible_SNPs, ref)

        for pos in sorted(all_results.keys()):
            r = all_results[pos][0]
            a = all_results[pos][1]
            results.append(f"{genome_basename}\t{chrom}\t{pos}\t{r}\t{a}\t")

        summary.append(
            f"{genome_basename}\t{chrom}\t{start_pos}\t{end_pos}\t{seq}\t{ref}\t"
            f"{';'.join(f'{pos}:{snp}' for pos, snp in SNPs.items()) if SNPs else '-'}\t"
            f"{len(SNPs) if SNPs else 0}\t"
            f"{';'.join(f'{pos}:{counts}' for pos, counts in variants.items()) if variants else '-'}\t"
            f"{len(variants) if variants else 0}\t"
            f"{';'.join(f'{pos}:{snp}' for pos, snp in possible_SNPs.items()) if possible_SNPs else '-'}\t"
            f"{false_seq_rate}\t"
        )

    genome_end_time = time.perf_counter()
    runtime = genome_end_time - genome_start_time
    print(f"Runtime for {genome_basename}: {runtime:.5f} seconds")

    # Clean up
    os.remove(f"/tmp/{genome_basename}.sam")


if output:
    output_dir = output
    os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist
    with open(f"{output_dir}/summary.tsv", "w") as o:
        o.write("\n".join(summary))

    with open(f"{output_dir}/results_tmp.tsv", "w") as o:
        o.write("\n".join(results))

    # Read the TSV file
    df = pd.read_csv(f"{output_dir}/results_tmp.tsv", sep="\t", header=0, index_col=False)

    # Transform the df
    pivot_df = df.pivot_table(
        values="Alternative", index=["Chromosome", "Position", "Reference"],
        columns="Genome", aggfunc='first'
    ).reset_index()

    # Rename columns
    pivot_df.columns = ["Chromosome", "Position", "Reference"] + [
        f"Alternative_{col}" for col in pivot_df.columns[3:]
    ]

    pivot_df.to_csv(f"{output_dir}/results.tsv", sep="\t", index=False)
    os.remove(f"{output_dir}/results_tmp.tsv")





