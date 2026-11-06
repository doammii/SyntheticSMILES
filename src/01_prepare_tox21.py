"""
Step 1 - prepare Tox21 and analyse structural alerts.

Outputs (into --outdir, default ./prepared):
  assays/<ASSAY>.csv          one file per assay: SMILES,label  (rows with a label only)
  alert_enrichment.csv        per assay: catalogue-level alert enrichment (odds ratio, p)
  alert_per_rule.csv          per assay x individual alert: enrichment, counts
  fragments/<ASSAY>.csv       concrete SMILES fragments to condition generation on
  summary.txt                 human-readable report

Nothing here needs a GPU; it takes a few minutes on CPU.

Why this step exists
--------------------
Known structural alerts (Brenk, NIH, PAINS) are external chemical knowledge: they
do not come from the Tox21 labels. If they are enriched among an assay's actives,
they carry label information that a generator can be conditioned on. If they are
not enriched, conditioning on them should not help - which makes those assays
negative controls. This script measures that enrichment so the downstream
experiment can be stratified by it.
"""

import argparse
import os
import urllib.request
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from scipy.stats import fisher_exact

RDLogger.DisableLog("rdApp.*")

TOX21_URL = ("https://raw.githubusercontent.com/GLambard/"
             "Molecules_Dataset_Collection/master/latest/tox21.csv")

CATALOGS = {
    "BRENK": FilterCatalogParams.FilterCatalogs.BRENK,
    "NIH": FilterCatalogParams.FilterCatalogs.NIH,
    "PAINS": FilterCatalogParams.FilterCatalogs.PAINS,
}


def load_tox21(cache):
    if not os.path.exists(cache):
        print(f"downloading Tox21 -> {cache}")
        urllib.request.urlretrieve(TOX21_URL, cache)
    df = pd.read_csv(cache)
    smi = [c for c in df.columns if c.lower() == "smiles"][0]
    assays = [c for c in df.columns
              if c not in (smi, "mol_id") and not c.lower().startswith("unnamed")]
    df = df.rename(columns={smi: "SMILES"})
    df["mol"] = df["SMILES"].map(lambda s: Chem.MolFromSmiles(s) if isinstance(s, str) else None)
    df = df[df["mol"].notna()].reset_index(drop=True)
    df["canon"] = df["mol"].map(Chem.MolToSmiles)
    return df, assays


def grow(mol, atom_idx, n_bonds):
    """Expand a set of atom indices outward by n_bonds."""
    s = set(atom_idx)
    for _ in range(n_bonds):
        s |= {nb.GetIdx() for i in list(s) for nb in mol.GetAtomWithIdx(i).GetNeighbors()}
    return s


def build_catalogs():
    out = {}
    for name, enum in CATALOGS.items():
        p = FilterCatalogParams()
        p.AddCatalog(enum)
        out[name] = FilterCatalog.FilterCatalog(p)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="prepared")
    ap.add_argument("--cache", default="tox21_original.csv")
    ap.add_argument("--min-rule-count", type=int, default=15,
                    help="minimum actives carrying an individual alert before testing it")
    ap.add_argument("--max-fragments", type=int, default=400,
                    help="max fragments written per assay")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "assays"), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "fragments"), exist_ok=True)

    df, assays = load_tox21(args.cache)
    print(f"parsed molecules: {len(df)} | assays: {len(assays)}")

    cats = build_catalogs()
    for name, cat in cats.items():
        df[name] = df["mol"].map(cat.HasMatch)
        print(f"  {name:6s} alert present in {100*df[name].mean():.1f}% of molecules")

    # ---- per-assay dataset files -------------------------------------------
    for a in assays:
        s = df[df[a].notna()][["SMILES", a]].rename(columns={a: "label"})
        s["label"] = s["label"].astype(int)
        s.to_csv(os.path.join(args.outdir, "assays", f"{a}.csv"), index=False)

    # ---- catalogue-level enrichment ----------------------------------------
    rows = []
    for a in assays:
        s = df[df[a].notna()]
        pos, neg = s[s[a] == 1], s[s[a] == 0]
        rec = {"assay": a, "n_labeled": len(s), "n_pos": len(pos),
               "pos_rate_pct": round(100 * len(pos) / len(s), 2),
               "imbalance_ratio": round(len(neg) / max(len(pos), 1), 1)}
        for name in cats:
            a11, a21 = int(pos[name].sum()), int(neg[name].sum())
            orr, p = fisher_exact([[a11, len(pos) - a11], [a21, len(neg) - a21]])
            rec[f"{name}_pos_pct"] = round(100 * a11 / len(pos), 1)
            rec[f"{name}_neg_pct"] = round(100 * a21 / len(neg), 1)
            rec[f"{name}_OR"] = round(float(orr), 3)
            rec[f"{name}_p"] = float(p)
        rows.append(rec)
    enr = pd.DataFrame(rows).sort_values("BRENK_OR", ascending=False)
    enr.to_csv(os.path.join(args.outdir, "alert_enrichment.csv"), index=False)

    # ---- individual-alert enrichment + concrete fragments -------------------
    # For every molecule, record which named alerts it carries and the concrete
    # SMILES of the matched atoms (this is what the generator is conditioned on).
    print("\nmatching individual alerts (this is the slow part)...")
    mol_rules = []          # list of dict: rule name -> set of fragment SMILES
    for m in df["mol"]:
        hits = {}
        for cat in cats.values():
            for e in cat.GetMatches(m):
                desc = e.GetDescription()
                try:
                    for fm in e.GetFilterMatches(m):
                        idx = [b for _, b in fm.atomPairs]
                        if len(idx) < 3:
                            continue
                        # the bare alert core, plus the core grown by one bond so
                        # the generator gets varied, molecule-specific context
                        for k in (0, 1):
                            s = grow(m, idx, k)
                            frag = Chem.MolFragmentToSmiles(m, atomsToUse=sorted(s),
                                                            canonical=True)
                            if frag:
                                hits.setdefault(desc, set()).add(frag)
                except Exception:
                    hits.setdefault(desc, set())
        mol_rules.append(hits)
    df["rules"] = mol_rules

    per_rule, frag_rows = [], defaultdict(list)
    for a in assays:
        s = df[df[a].notna()]
        pos, neg = s[s[a] == 1], s[s[a] == 0]
        pos_counts = Counter(r for h in pos["rules"] for r in h)
        neg_counts = Counter(r for h in neg["rules"] for r in h)
        cand = []
        for rule, k in pos_counts.items():
            if k < args.min_rule_count:
                continue
            kn = neg_counts.get(rule, 0)
            orr, p = fisher_exact([[k, len(pos) - k], [kn, len(neg) - kn]])
            cand.append((rule, k, kn, float(orr), float(p)))
        # Benjamini-Hochberg over the rules tested for this assay
        cand.sort(key=lambda r: r[4])
        mt = len(cand)
        for i, (rule, k, kn, orr, p) in enumerate(cand, start=1):
            q = min(1.0, p * mt / i)
            per_rule.append({"assay": a, "alert": rule, "n_pos_with": k,
                             "n_neg_with": kn,
                             "pos_pct": round(100 * k / len(pos), 1),
                             "neg_pct": round(100 * kn / len(neg), 1),
                             "odds_ratio": round(orr, 3), "p": p, "q_BH": q,
                             "enriched": bool(orr > 1 and q < 0.05)})
        # fragments: from actives carrying a significantly enriched alert
        keep = {r["alert"] for r in per_rule
                if r["assay"] == a and r["enriched"]}
        seen = set()
        for h in pos["rules"]:
            for rule, frags in h.items():
                if rule not in keep:
                    continue
                for f in frags:
                    if f not in seen and Chem.MolFromSmiles(f) is not None:
                        seen.add(f)
                        frag_rows[a].append({"fragment": f, "alert": rule})
        pd.DataFrame(frag_rows[a][: args.max_fragments]).to_csv(
            os.path.join(args.outdir, "fragments", f"{a}.csv"), index=False)

    pd.DataFrame(per_rule).to_csv(os.path.join(args.outdir, "alert_per_rule.csv"), index=False)

    # ---- report -------------------------------------------------------------
    lines = ["Tox21 preparation summary", "=" * 70, ""]
    lines.append(f"molecules parsed : {len(df)}")
    lines.append(f"assays           : {len(assays)}")
    lines.append("")
    lines.append("Per-assay imbalance and BRENK alert enrichment")
    lines.append(f"{'assay':16s} {'n':>6s} {'pos':>6s} {'pos%':>6s} {'1:x':>7s} "
                 f"{'OR':>6s} {'p':>10s} {'frags':>6s}")
    for _, r in enr.iterrows():
        nf = len(frag_rows[r["assay"]])
        lines.append(f"{r['assay']:16s} {int(r['n_labeled']):6d} {int(r['n_pos']):6d} "
                     f"{r['pos_rate_pct']:6.1f} {r['imbalance_ratio']:7.1f} "
                     f"{r['BRENK_OR']:6.2f} {r['BRENK_p']:10.2e} {nf:6d}")
    lines += ["", "Reading this table:",
              "  OR > 1 and small p  -> alerts carry label signal; alert-conditioned",
              "                         generation is expected to help here.",
              "  OR ~ 1              -> neutral assay.",
              "  OR < 1              -> alerts are DEPLETED among actives; these assays",
              "                         are negative controls - alert conditioning should",
              "                         NOT help, and if it appears to, something is wrong.", ""]
    txt = "\n".join(lines)
    open(os.path.join(args.outdir, "summary.txt"), "w").write(txt)
    print("\n" + txt)
    print(f"written to {args.outdir}/")


if __name__ == "__main__":
    main()
