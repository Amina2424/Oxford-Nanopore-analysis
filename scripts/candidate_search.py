#!/usr/bin/env python3

import pandas as pd

def find_regions(bedmethyl_file, output_file, min_methylation=43, max_methylation=67, min_cpg=5, max_length=6000):
    df = pd.read_csv(bedmethyl_file, sep="\t", header=None)
    df.columns = [
    "chrom", "start", "end", "mod_code", "score", "strand",
    "start_compat", "end_compat", "color", "Nvalid_cov",
    "methylation_level", "Nmod", "Ncanonical", "Nother_mod",
    "Ndelete", "Nfail", "Ndiff", "N_local"]

    regions = []
    current_region = []
    current_chrom = None


    for index, row in df.iterrows():
        if current_chrom != row['chrom']:
            if current_region:
                regions.append(current_region)
            current_region = []
            current_chrom = row['chrom']


        if min_methylation <= row['methylation_level'] <= max_methylation:
            current_region.append(row)
        else:
            if len(current_region) >= min_cpg:
                regions.append(current_region)
            current_region = []


    if len(current_region) >= min_cpg:
        regions.append(current_region)


    filtered_regions = []
    for region in regions:
        start = region[0]['start']
        end = region[-1]['end']
        if end - start <= max_length:
            filtered_regions.append(region)


    with open(output_file, 'w') as f:
        f.write("chrom\tstart\tend\tcomputed_methylation_level\tall_methylation_levels\n")
        for region in filtered_regions:
            chrom = region[0]['chrom']
            start = region[0]['start']
            end = region[-1]['end']
            length = region[-1]['end'] - region[0]['start'] + 1
            methylation_levels = [r['methylation_level'] for r in region]
            computed_methylation_level = sum(methylation_levels) / len(methylation_levels)
            f.write(f"{chrom}\t{start}\t{end}\t{length}\t{computed_methylation_level}\t{','.join(map(str, methylation_levels))}\n")

bedmethyl_file = input()
output_file = 'predicted_monoallelic_regions.bed'
find_regions(bedmethyl_file, output_file)
