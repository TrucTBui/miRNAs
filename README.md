# miRNAs

## bam2seq.py

This script identifies variants within specified regions of a provided genome file (BAM format).

Features
- Extracts reads from the BAM file for a given location.
- Trims reads to match the specified region.
- Identifies variants within the trimmed reads and their frequencies.
- Outputs results to the console or a user-defined file.
- Measures and reports runtime for each processed genome file.

Usage: \
python3 bam2seq.py -g <genome_file1> <genome_file2> ... -l <location1> <location2> ... [-o <output_file>]

Arguments: \
-g: Required. Path to one or more genome files in BAM format. Provide separate paths for multiple files. \
-l: Required. Location(s) within the genome to analyze. Use format <chromosome>:<start_position>-<end_position> for each location. Provide separate locations with spaces. \
-o (Optional): Path to an output file to store the results. If not provided, results will be printed to the console. 

Example:\
python3 bam2seq.py -g human_genome.bam -l chr1:1000-1005 chr2:2000-2010 -o variants.txt \
This command will analyze regions chr1:1000-1005 and chr2:2000-2010 in the human_genome.bam file and write the results (sequence and variant information) to the variants.txt file. 

NOTE: Please ensure that the chromosome names in your location specifications match the format used in your BAM files. Some BAM files use "chr1", "chr2", etc., while others use "1", "2", etc. 