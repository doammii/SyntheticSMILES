"""
Step 2 - the experiment.

Runs one Tox21 assay under two evaluation protocols and several augmentation
conditions, with repeated seeds, and writes per-seed and summary metrics.

Protocols
---------
clean         split the ORIGINAL data first, then augment the TRAINING split only.
              The test set contains no derived or generated molecule. This is the
              correct protocol.
contaminated  augment the whole dataset first, then split. Derived and generated
              molecules land in validation and test too. This reproduces the
              common error and is included deliberately so the inflation can be
              measured, not as a mistake.

Conditions (all add the same number of rows, so they differ only in source)
--------------------------------------------------------------------------
baseline   no augmentation
dup        duplication of real minority molecules
randsmi    randomized (non-canonical) SMILES of real minority molecules
llamol     generated molecules from --synthetic

With --synthetic, a `parent` column (the molecule the generation was conditioned
on) is used under the clean protocol to keep only generations whose parent is in
the training split. Without it, all generated molecules are used and the script
says so - that is a weaker guarantee.

Example
-------
  python 02_run_experiments.py --data prepared/assays/SR-ARE.csv \
      --outdir out/SR-ARE --protocol clean --seeds 0 1 2 3 4
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from sklearn.metrics import (
    average_precision_score, confusion_matrix, matthews_corrcoef, roc_auc_score,
)
from sklearn.model_selection import train_test_split

RDLogger.DisableLog("rdApp.*")

MODEL_NAME = "DeepChem/ChemBERTa-10M-MLM"
LR = 1e-5
EPOCHS = 50
PATIENCE = 10
BATCH_SIZE = 32
MAX_SEQ_LEN = 128


def canon(s):
    if not isinstance(s, str) or not s.strip():
        return None
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else None


def read_table(path, need_parent=False):
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    smi = cols.get("smiles") or list(df.columns)[0]
    lab = cols.get("label") or cols.get("labels")
    out = df.rename(columns={smi: "SMILES", lab: "label"})
    keep = ["SMILES", "label"] + (["parent"] if need_parent and "parent" in df.columns else [])
    out = out[keep].copy()
    out["label"] = pd.to_numeric(out["label"], errors="coerce")
    out = out.dropna(subset=["label"])
    out["label"] = out["label"].astype(int)
    out["canon"] = out["SMILES"].map(canon)
    out = out.dropna(subset=["canon"])
    if "parent" in out.columns:
        out["parent_canon"] = out["parent"].map(canon)
    return out.reset_index(drop=True)


def dedupe(df, what):
    bad = df.groupby("canon")["label"].nunique()
    bad = set(bad[bad > 1].index)
    if bad:
        print(f"[{what}] dropping {len(bad)} molecules with conflicting labels")
        df = df[~df["canon"].isin(bad)]
    n0 = len(df)
    df = df.drop_duplicates(subset="canon").reset_index(drop=True)
    if n0 != len(df):
        print(f"[{what}] removed {n0 - len(df)} duplicate molecules")
    return df


def duplicate_minority(pool, n_add, minority, rng):
    sub = pool[pool["label"] == minority]
    idx = rng.choice(len(sub), size=n_add, replace=True)
    return sub.iloc[idx][["SMILES", "label"]].copy()


def randomized_smiles(pool, n_add, minority, rng):
    sub = pool[pool["label"] == minority]
    out, guard = [], 0
    while len(out) < n_add and guard < n_add * 50:
        guard += 1
        m = Chem.MolFromSmiles(sub.iloc[rng.integers(len(sub))]["SMILES"])
        if m is None:
            continue
        out.append({"SMILES": Chem.MolToSmiles(m, canonical=False, doRandom=True),
                    "label": minority})
    return pd.DataFrame(out)


def _cuda():
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


def run_one(train, valid, test, seed, tag, outdir):
    from simpletransformers.classification import ClassificationModel

    args = {
        "num_train_epochs": EPOCHS, "train_batch_size": BATCH_SIZE,
        "eval_batch_size": BATCH_SIZE, "learning_rate": LR,
        "max_seq_length": MAX_SEQ_LEN, "manual_seed": seed,
        "use_early_stopping": True, "early_stopping_patience": PATIENCE,
        "early_stopping_metric": "eval_loss", "early_stopping_metric_minimize": True,
        "evaluate_during_training": True, "evaluate_each_epoch": True,
        "save_best_model": True, "save_model_every_epoch": False,
        "save_eval_checkpoints": False, "save_steps": -1,
        "overwrite_output_dir": True, "auto_weights": True, "silent": True,
        "output_dir": os.path.join(outdir, "_work", f"{tag}_s{seed}"),
        "best_model_dir": os.path.join(outdir, "_work", f"{tag}_s{seed}", "best"),
        "use_multiprocessing": False, "use_multiprocessing_for_evaluation": False,
    }

    def fmt(d):
        return d[["SMILES", "label"]].rename(columns={"SMILES": "text", "label": "labels"})

    model = ClassificationModel("roberta", MODEL_NAME, args=args, use_cuda=_cuda())
    model.train_model(fmt(train), eval_df=fmt(valid))
    if os.path.isdir(args["best_model_dir"]):
        model = ClassificationModel("roberta", args["best_model_dir"], args=args,
                                    use_cuda=_cuda())

    def probs(d):
        _, raw = model.predict(list(d["SMILES"]))
        raw = np.asarray(raw, dtype=float)
        e = np.exp(raw - raw.max(axis=1, keepdims=True))
        return (e / e.sum(axis=1, keepdims=True))[:, 1]

    pv, pt = probs(valid), probs(test)
    thr = best_threshold(valid["label"].values, pv)
    return {"at_0.5": metrics(test["label"].values, pt, 0.5),
            "at_valid_threshold": metrics(test["label"].values, pt, thr),
            "threshold": float(thr)}


def best_threshold(y, p):
    g = np.linspace(0.05, 0.95, 91)
    return float(g[int(np.argmax([matthews_corrcoef(y, (p >= t).astype(int)) for t in g]))])


def metrics(y, p, thr):
    pred = (p >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
            "sensitivity": float(sens), "specificity": float(spec),
            "balanced_accuracy": float((sens + spec) / 2),
            "mcc": float(matthews_corrcoef(y, pred)),
            "auroc": float(roc_auc_score(y, p)),
            "auprc": float(average_precision_score(y, p))}


def main():
    global EPOCHS, PATIENCE, BATCH_SIZE, MAX_SEQ_LEN
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="assay CSV with SMILES,label")
    ap.add_argument("--synthetic", help="generated molecules: SMILES,label[,parent]")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--protocol", choices=["clean", "contaminated"], default="clean")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--conditions", nargs="+",
                    default=["baseline", "dup", "randsmi", "llamol"])
    ap.add_argument("--n-add", type=int, default=0,
                    help="rows to add per condition; 0 = match the generated set, "
                         "or 2x the minority count when there is none")
    ap.add_argument("--test-size", type=float, default=0.2)
    ap.add_argument("--valid-size", type=float, default=0.2)
    ap.add_argument("--epochs", type=int, default=EPOCHS)
    ap.add_argument("--patience", type=int, default=PATIENCE)
    ap.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    ap.add_argument("--max-seq-len", type=int, default=MAX_SEQ_LEN)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    EPOCHS, PATIENCE = a.epochs, a.patience
    BATCH_SIZE, MAX_SEQ_LEN = a.batch_size, a.max_seq_len

    original = dedupe(read_table(a.data), "data")
    minority = int(original["label"].value_counts().idxmin())
    n_min = int((original["label"] == minority).sum())

    synthetic = None
    if a.synthetic and os.path.exists(a.synthetic):
        synthetic = dedupe(read_table(a.synthetic, need_parent=True), "synthetic")
        synthetic = synthetic[~synthetic["canon"].isin(set(original["canon"]))].reset_index(drop=True)
        if "parent_canon" not in synthetic.columns:
            print("WARNING: no `parent` column in the generated file. Under the clean "
                  "protocol every generated molecule will be used, so generations "
                  "derived from test molecules cannot be excluded.")
    if "llamol" in a.conditions and synthetic is None:
        a.conditions = [c for c in a.conditions if c != "llamol"]
        print("no --synthetic given; skipping the llamol condition")

    n_add = a.n_add or (len(synthetic) if synthetic is not None else 2 * n_min)

    print(f"device    : {'cuda' if _cuda() else 'cpu'} | protocol {a.protocol}")
    print(f"assay data: {len(original)} {dict(original['label'].value_counts())} "
          f"(minority class = {minority}, n = {n_min})")
    if synthetic is not None:
        print(f"generated : {len(synthetic)} {dict(synthetic['label'].value_counts())}")
    print(f"rows added per augmented condition: {n_add}")

    rows = []
    for seed in a.seeds:
        rng = np.random.default_rng(seed)

        if a.protocol == "clean":
            trval, test = train_test_split(original, test_size=a.test_size,
                                           stratify=original["label"], random_state=seed)
            train, valid = train_test_split(trval, test_size=a.valid_size,
                                            stratify=trval["label"], random_state=seed)
            pool = train                      # augmentation draws from the train split
            if synthetic is not None:
                if "parent_canon" in synthetic.columns:
                    keep = synthetic["parent_canon"].isin(set(train["canon"]))
                    syn_use = synthetic[keep]
                else:
                    syn_use = synthetic
            else:
                syn_use = None
        else:
            # contaminated: build the augmented table first, then split it
            pool = original
            syn_use = synthetic

        if seed == a.seeds[0]:
            print(f"\nsplit check (seed {seed}): "
                  f"{'train %d / valid %d / test %d' % (len(train), len(valid), len(test)) if a.protocol=='clean' else 'split happens after augmentation'}")

        for tag in a.conditions:
            if tag == "baseline":
                add = None
            elif tag == "dup":
                add = duplicate_minority(pool, n_add, minority, rng)
            elif tag == "randsmi":
                add = randomized_smiles(pool, n_add, minority, rng)
            elif tag == "llamol":
                if syn_use is None or len(syn_use) == 0:
                    print(f"seed {seed}: no usable generated molecules, skipping llamol")
                    continue
                take = min(n_add, len(syn_use))
                if take < n_add:
                    print(f"  WARNING: seed {seed} has only {len(syn_use)} usable "
                          f"generated molecules but {n_add} rows are needed, so the "
                          f"llamol condition is NOT volume-matched to dup/randsmi "
                          f"here. Generate more molecules and re-run before "
                          f"comparing.", flush=True)
                add = syn_use.sample(n=take, random_state=seed)[["SMILES", "label"]]
            else:
                raise SystemExit(f"unknown condition {tag}")

            if a.protocol == "clean":
                tr = train[["SMILES", "label"]]
                if add is not None:
                    tr = pd.concat([tr, add], ignore_index=True)
                va, te = valid, test
            else:
                full = original[["SMILES", "label"]]
                if add is not None:
                    full = pd.concat([full, add], ignore_index=True)
                trval_c, te = train_test_split(full, test_size=a.test_size,
                                               stratify=full["label"], random_state=seed)
                tr, va = train_test_split(trval_c, test_size=a.valid_size,
                                          stratify=trval_c["label"], random_state=seed)

            leak = 0
            if add is not None and a.protocol == "contaminated":
                leak = int(te["SMILES"].isin(set(add["SMILES"])).sum())
            print(f"\n--- {a.protocol} | seed {seed} | {tag} | train {len(tr)} "
                  f"(minority {int((tr['label']==minority).sum())}) | test {len(te)} "
                  f"(minority {int((te['label']==minority).sum())}"
                  f"{', of which added: %d' % leak if leak else ''}) ---", flush=True)

            res = run_one(tr, va, te, seed, f"{a.protocol}_{tag}", a.outdir)
            for mode in ("at_0.5", "at_valid_threshold"):
                rows.append({"assay": os.path.basename(a.data).replace(".csv", ""),
                             "protocol": a.protocol, "seed": seed, "condition": tag,
                             "mode": mode, "threshold": res["threshold"],
                             "added_rows_in_test": leak, **res[mode]})
            print(json.dumps(res["at_valid_threshold"]), flush=True)
            pd.DataFrame(rows).to_csv(os.path.join(a.outdir, "results_raw.csv"), index=False)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.outdir, "results_raw.csv"), index=False)
    cols = ["sensitivity", "specificity", "balanced_accuracy", "mcc", "auroc", "auprc"]
    s = df.groupby(["mode", "condition"])[cols].agg(["mean", "std"]).round(4)
    s.to_csv(os.path.join(a.outdir, "results_summary.csv"))
    print("\n============ SUMMARY (mean / sd over seeds) ============")
    print(s.to_string())
    print(f"\nwritten to {a.outdir}/")


if __name__ == "__main__":
    main()
