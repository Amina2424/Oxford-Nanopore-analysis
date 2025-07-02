#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import bootstrap

path = input('файл: ')
out = input('путь выхода:')
df = pd.read_csv(path, sep="\t", header=None)
df.columns = [
    "chrom", "start", "end", "mod_code", "score", "strand",
    "start_compat", "end_compat", "color", "Nvalid_cov",
    "percent_modified", "Nmod", "Ncanonical", "Nother_mod",
    "Ndelete", "Nfail", "Ndiff", "N_local"
]

def bootstrap_ci(data, n_resamples=100000, confidence_level=0.99):
    res = bootstrap((data,), np.mean, n_resamples=n_resamples, confidence_level=confidence_level, random_state=42)
    mean = round(np.mean(res.bootstrap_distribution), 2)

    return round(res.confidence_interval.low,3), round(res.confidence_interval.high, 3), mean

region_data = df.query('mod_code == "m"')
ci_low, ci_high, mean_v = bootstrap_ci(region_data["percent_modified"].values)

pd.DataFrame([[ci_low, ci_high, mean_v]], columns=["CI_lower", "CI_upper", 'Mean']).to_csv(out, sep='\t', header=True, index=False)

print(f"Тук-тук, готово")
