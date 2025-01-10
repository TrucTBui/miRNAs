import subprocess
import argparse
import os


input_genome_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/genomes"

def run(gene_name, ID):
    for filename in os.listdir(input_genome_folder):
        file_path = os.path.join(input_genome_folder, filename)
        output_folder = os.path.join(f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/{gene_name}",
                                     os.path.splitext(os.path.basename(filename))[0])
        cmd = f"python3 /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/bam_fa2seq.py --genome {file_path} --reference /mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz --location /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/{gene_name}/annot_{ID}.tsv --output {output_folder}/"
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)

gene_list = {"ENSG00000139618":"BRCA2", "ENSG00000160791":"CCR5","ENSG00000240972":"MIF"}
#gene_list = {"ENSG00000139618":"BRCA2"}
#gene_list ={"ENSG00000160791":"CCR5"}
#gene_list ={"ENSG00000240972":"MIF"}

for id, gene_name in gene_list.items():
    run(gene_name, id)