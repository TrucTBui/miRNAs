# miRNAs

## bam_fa2seq.py

This script analyzes sequencing data from multiple genomes (provided as BAM files) to identify variants within specified regions of interest. It compares these variants to a reference genome provided in FASTA format.

**Features**:

- Extracts reads from BAM files for specific genomic regions.
- Identifies variants (SNPs) within the reads compared to the reference genome.
- Calculates an approximate false sequencing error rate.
- Generates two output files (configurable):
  - summary.tsv: Summarizes variants and SNPs found for each region across all genomes.
  - results.tsv: Detailed information about each variant call, including reference and alternative alleles for each genome analyzed.

**Requirements**:
- Python 3+
- pandas library
- samtools 

**Usage**:

python3 bam_fa2seq.py -g <genomes.txt> -r <reference.fa> -l <location.txt> -o <output_dir> 
Arguments:

-g: Path to a text file containing paths to BAM files (one per line).
-r: Path to the reference genome FASTA file.
-l: Path to a text file containing genomic regions to analyze (one region per line, format: chr:start-end).
-o (Optional): Path to the output directory. Defaults to the current working directory.

Requirements for input Files:
- genomes.txt: A text file containing paths to each BAM file, one per line.
- location.txt: A text file containing genomic regions to analyze, one per line. Each line should be formatted as chromosome:start-end (e.g., chr1:1000-2000).

Output Files:
- summary.tsv: A tab-delimited file summarizing variants and SNPs found for each region across all genomes.
- results.tsv: A tab-delimited file containing detailed information about each variant call, including reference and alternative alleles for each genome analyzed. (Written to a temporary file first, then transformed and written to the final output file.)

Notes:
- Please ensure that the chromosome names in your location specifications match the format used in your BAM files. Some BAM files use "chr1", "chr2", etc., while others use "1", "2", etc. 
- This script uses temporary files during execution, which are automatically cleaned up.
- The script calculates an approximate false sequencing error rate based on the identified variants.
- This script provides a basic framework for variant calling and analysis. You can modify it to suit your specific needs, such as filtering variants based on quality scores or incorporating additional analysis steps.
- This script relies on samtools, which is installed and accessible from the LMU Bioinformatics chair's internal environment. You may need to adjust the samtools command if used in a different setup.