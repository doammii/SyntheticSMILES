"""Combine all results_raw.csv files and print mean +/- sample std (ddof=1) per assay, protocol, condition."""
import glob, os, pandas as pd
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
files = glob.glob(os.path.join(root, "results", "raw", "*", "results_raw.csv"))
df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
df.to_csv(os.path.join(root, "results", "all_results_raw.csv"), index=False)
m = ["sensitivity", "specificity", "balanced_accuracy", "mcc", "auroc", "auprc"]
g = df.groupby(["assay", "protocol", "condition", "mode"])[m].agg(["mean", "std"])  # std = ddof=1
g.to_csv(os.path.join(root, "results", "summary_mean_std.csv"))
print(g.round(3).to_string())
print("rows:", len(df), "| files:", len(files))
