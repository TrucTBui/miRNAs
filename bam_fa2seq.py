"""
version 06.12.2024
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
    variants = {}
    all_results = {}
    processed_results = {}
    SNPs = {}
    haploids = {}
    counter = 0  # for iterating the seq from ref genome

    for pos in sorted(variant_dict.keys()):
        ref = reference[counter]
        all_results[pos] = [ref]  # 0.index contains ref, 1.index contains alternative

        bases = variant_dict[pos]
        frequency = Counter(bases)  # Count how frequent each variant is
        unique_bases = set(bases)
        all_results[pos].append("/".join(unique_bases)) # 0.index contains ref, 1.index contains alternative

        # Find variances in reads
        if len(unique_bases) > 1:  # More than one unique nucleotide means a variant
            sorted_bases = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
            ordered_bases = [base for base, count in sorted_bases]

            variants[pos] = dict(sorted_bases)

            most_frequent_base = ordered_bases[0]  # Extract the most significant base
            ratio = (frequency[most_frequent_base] / sum(frequency.values()))

            # Process the reads results: The most frequent base is supposed to be the correct one
            if ratio >= 0.75 or (len(unique_bases) > 2 and ratio >= 0.5):  # 3 or 4 bases detected when reading
                processed_results[pos] = most_frequent_base
                seq += most_frequent_base

                # compare the most frequent base read at position with the ref
                if not most_frequent_base == ref:  # SNP found
                    SNPs[pos] = f"{ref}/{most_frequent_base}"
            else:
                processed_results[pos] = "/".join(unique_bases)
                seq += "(" + "/".join(ordered_bases) + ")"

                haploids[pos] = "/".join(ordered_bases)


        else:  # Only one base at that position
            base = bases[0]
            seq += base
            processed_results[pos] = base

            if not base == ref:  # SNP found
                SNPs[pos] = f"{ref}/{base}"

        counter += 1

    return all_results, processed_results, seq, variants, SNPs, haploids

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

def transform_location_input(path):
    locations= []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.split()
                chromosome = fields[3]
                if chromosome.startswith("chr"):
                    chromosome = chromosome[3:]
                start = fields[5]
                end = fields[6]
                locations.append(f"{chromosome}:{start}-{end}")
    with open(f"{os.path.dirname(path)}/locations_tranformed.txt", 'w') as w:
         for location in locations:
             w.write(f"{location}\n")
    return locations

def calc_false_sequencing_rate(variants:dict,processed_results, ref:str):
    if len(variants) > 0:
        error_num = 0
        for pos in variants.keys():
            if "/" not in processed_results[pos]:
                error_num+=1
        return round(error_num/len(ref), 5)
    else:
        return 0.0

def print_output(output_path):
    os.makedirs(output_path, exist_ok=True)  # Create the directory if it doesn't exist
    with open(f"{output_path}/summary.tsv", "w") as o:
        o.write("\n".join(summary))

    with open(f"{output_path}/results_tmp.tsv", "w") as o:
        o.write("\n".join(results))

    # Read the TSV file
    df = pd.read_csv(f"{output_path}/results_tmp.tsv", sep="\t", header=0, index_col=False)

    # Transform the df
    pivot_df = df.pivot_table(
        values=["Alternative", "Variants_Freq", "Processed_Alt"], index=["Chromosome", "Position", "Reference"],
        columns="Genome", aggfunc='first'
    ).reset_index()

    pivot_df.columns = ['Chromosome', 'Position', 'Reference'] + [
        f"{col[0]}_{col[1]}" for col in pivot_df.columns[3:]
    ]


    # Rearrange the columns
    columns = ['Chromosome', 'Position', 'Reference']
    other_columns = [col for col in pivot_df.columns if col not in columns]
    genome_names = set(a.split("_")[1] for a in other_columns)
    for name in genome_names:
        columns.extend([col for col in pivot_df.columns if col.endswith(name)])
    pivot_df = pivot_df[columns]

    pivot_df['Conclusion'] = 'Normal'

    problematic_pos = []
    # raw_base_cols = [col for col in columns if col.startswith("Alternative")]
    final_base_cols = [col for col in columns if col.startswith("Processed_Alt")]

    for index, row in pivot_df.iterrows():
        ref = row.Reference
        alts = []
        """
        for col in raw_base_cols:
            if row[col].count('/') == 3:  # Problematic: All 4 bases were read at the position
                pivot_df.loc[index, 'Conclusion'] = 'Problem'
                problematic_pos.append(f"{row.Chromosome}:{row.Position}")
                break
        """
        for col in final_base_cols:
            alt = row[col]
            alts.append(alt)
        if len(set(alts)) > 1 or any('/' in alt for alt in alts):  # Conflicts in read
            if all('/' in alt for alt in alts) or all('/' not in alt for alt in alts):
                pivot_df.loc[index, 'Conclusion'] = 'Haploid'
            else:  # Problematic: Haploid in one genome and no-halpoid in the other one
                pivot_df.loc[index, 'Conclusion'] = 'Problem'
                problematic_pos.append(f"{row.Chromosome}:{row.Position}")
        elif ref != alts[0]:  # Base is clear and different from the ref
            pivot_df.loc[index, 'Conclusion'] = 'SNP'


    pivot_df.to_csv(f"{output_path}/results.tsv", sep="\t", index=False)

    with open(f"{output_path}/problems.csv", "w") as o:
        if len(problematic_pos) > 0:
            grouped_regions = group_problematic_positions(problematic_pos)
            o.write("\n".join(grouped_regions))
        else:
            o.write("No problematic positions detected.")

    os.remove(f"{output_path}/results_tmp.tsv")


def group_problematic_positions(problematic_pos, distance_threshold=5):
    grouped_regions = []
    current_region = None

    for pos in sorted(problematic_pos):
        chrom, position = pos.split(":")
        position = int(position)

        if current_region is None:
            current_region = [chrom, position, position]  # Start a new region
        elif (chrom == current_region[0]) and (position <= current_region[2] + distance_threshold):
            # Extend the current region if within the threshold
            current_region[2] = position
        else:
            # Save the completed region and start a new one
            grouped_regions.append([current_region[0], current_region[1], current_region[2]])
            current_region = [chrom, position, position]

    # Add the last region
    if current_region is not None:
        grouped_regions.append([current_region[0], current_region[1], current_region[2]])

    # Merging
    def merge_regions(regions):
        merged = []
        # regions = sorted(regions, key=lambda x: x[1])  # Sort by start position

        while regions:
            current = regions.pop(0)
            # Try to merge with the next region if close enough
            # Check if the end postition of the current region is close to start position of the next region
            while regions and (regions[0][0] == current[0]) and (regions[0][1] <= current[2] + distance_threshold):
                current[2] = max(current[2], regions[0][2])  # Extend the current region
                regions.pop(0)  # Remove the merged region

            merged.append(f"{current[0]}:{current[1]}-{current[2]}")

        return merged

    # Perform the merge
    return merge_regions(grouped_regions)


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

"""
# The input file is already a csv file, each line containing the position
if os.path.isfile(location_path):  # Input is a file
    with open(location_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                locations.append(line.strip())
else:
    raise Exception("Location file not found or not correct")

ref_genome = extract_ref_genome(ref_path, location_path)
"""

if os.path.isfile(location_path):  # Input is the annotation file from ISAR
    locations = transform_location_input(location_path)
    transformed_path = f"{os.path.dirname(location_path)}/locations_tranformed.txt"
    ref_genome = extract_ref_genome(ref_path, transformed_path)

else:
    raise Exception("Location file not found or not correct")


results = [f"Genome\tChromosome\tPosition\tReference\tAlternative\tVariants_Freq\tProcessed_Alt"]

summary = [f"Genome\tChromosome\tStart\tEnd\tProcessed_Sequence\tReference\tSNP\tSNP_Freq\tSequencing_Variance_Freq\tHaploid\tApprox_sequencing_rate"]

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

        all_results, processed_results, seq, variants, SNPs, haploids = find_variants_per_pos(variant_dict, ref)

        chrom, pos_range = location.split(":")
        start_pos, end_pos = map(int, pos_range.split("-"))

        false_seq_rate = calc_false_sequencing_rate(variants,processed_results, ref)

        for pos in sorted(all_results.keys()):
            r = all_results[pos][0]
            a = all_results[pos][1]
            p = processed_results[pos]  # Final base
            if pos in variants.keys():
                freq = ','.join(f"{key}:{value}" for key, value in variants[pos].items())
            else:
                freq = "-"
            results.append(f"{genome_basename}\t{chrom}\t{pos}\t{r}\t{a}\t{freq}\t{p}")

        summary.append(
            f"{genome_basename}\t{chrom}\t{start_pos}\t{end_pos}\t{seq}\t{ref}\t"
            f"{';'.join(f'{pos}:{snp}' for pos, snp in SNPs.items()) if SNPs else '-'}\t"
            f"{len(SNPs) if SNPs else 0}\t"
            f"{len(variants) if variants else 0}\t"
            f"{';'.join(f'{pos}:{snp}' for pos, snp in haploids.items()) if haploids else '-'}\t"
            f"{false_seq_rate}\t"
        )

    genome_end_time = time.perf_counter()
    runtime = genome_end_time - genome_start_time
    print(f"Runtime for {genome_basename}: {runtime:.5f} seconds")

    # Clean up
    os.remove(f"/tmp/{genome_basename}.sam")


print_output(output)







