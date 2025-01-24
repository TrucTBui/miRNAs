import pandas as pd
def find_diff_in_results_file(result, output):
    result_df = pd.read_csv(result, sep="\t", header=0, index_col=False,low_memory=False)
    differences=[]
    final_base_cols = [col for col in result_df.columns if col.startswith("Processed_Alt")]

    for index, row in result_df.iterrows():
        ref = row.Reference  # Reference base

        alts = [row[col] for col in final_base_cols]  # Collect all bases from the family


        if len(set(alts)) > 1 or any('/' in alt for alt in alts):  # Diff in read
           differences.append(row)

    diff_df = pd.DataFrame(differences)

    diff_df.to_csv(output, sep="\t", index=False)


find_diff_in_results_file("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/CCR5/results.tsv",
                          "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/CCR5/results_filtered.tsv")