"""
Step 3 - generate molecules with LLaMOL, conditioned on structural alerts.

Run this from inside the llamol repository directory, because it imports the
repo's own `sample` module and loads `out/llama2-M-Full-RSS-Canonical.pt`:

    cd /path/to/llamol-main
    python /path/to/03_generate_llamol.py --smoke-test
    python /path/to/03_generate_llamol.py \
        --assay-csv .../prepared/assays/SR-ARE.csv \
        --fragments .../prepared/fragments/SR-ARE.csv \
        --mode alert --out .../generated/SR-ARE_alert.csv

Three generation modes
----------------------
  whole   context = the whole minority molecule. Reproduces what the earlier
          experiments did, as a comparison arm. The only mode available for an
          assay whose alert set is empty (NR-AR, NR-AR-LBD, NR-ER-LBD).
  alert   context = a structural-alert fragment found in minority molecules.
          External chemical knowledge, independent of the labels.
  mixed   half of each.

What is different from the earlier pipeline
-------------------------------------------
1. The property conditions (logP, SA score, molecular weight / 100) are computed
   from the parent molecule instead of being left at the demo defaults
   (0.0 / 2.0 / 3.0), which were used unchanged for every molecule and every
   dataset before.
2. Every generated molecule records its parent and its conditioning fragment, so
   downstream splits can keep only generations whose parent is in the training
   set. This makes evaluation leakage impossible by construction.
3. Generated molecules are filtered: RDKit-valid, deduplicated, absent from the
   assay data, and - in alert mode - required to retain the conditioning
   substructure.

Revision 2026-10-07
-------------------
Three things were wrong in the first version and are fixed here.

a) A fragment that matches no active molecule used to fall back to using the
   fragment string itself as the parent. Such a parent is never in the training
   split, so step 2 discarded every generation derived from it under the clean
   protocol - silently. These fragments are now skipped and counted. They are
   an artefact of charge notation in the step-1 alert patterns (for example
   `O=[N+][O-]`, `C=CC(=O)[O]`), which do not substructure-match the molecules
   they were extracted from. 21 of 147 fragments in SR-ARE and 17 of 105 in
   SR-MMP are affected.

b) Only the first matching active was used as the parent for a fragment, so all
   generations from that fragment shared one parent. Whenever that parent landed
   in validation or test, the whole block vanished under the clean protocol, and
   the usable count swung wildly between seeds. Each fragment is now paired with
   up to --parents-per-fragment distinct parents, drawn in shuffled order, which
   both spreads the parents across the splits and multiplies the number of
   contexts for assays with few fragments.

c) --max-contexts silently truncated the context list at 400. It now defaults to
   no limit; pass a positive number to cap it.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors

RDLogger.DisableLog("rdApp.*")

try:
    from rdkit.Chem import RDConfig
    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer
except Exception:
    sascorer = None


def props(mol):
    """LLaMOL's three context values. mol_weight is Daltons / 100."""
    logp = float(Crippen.MolLogP(mol))
    sa = float(sascorer.calculateScore(mol)) if sascorer else 2.5
    mw = float(Descriptors.MolWt(mol)) / 100.0
    return logp, sa, mw


def make_sampler(ckpt, seed, device=None):
    import torch
    from sample import Sampler
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = "float16" if device == "cuda" else "float32"
    return Sampler(load_path=ckpt, device=device, seed=seed, dtype=dtype,
                   compile=False), device


def generate_for(sampler, device, context_smi, logp, sa, mw, n, temperature, top_k,
                 max_new_tokens):
    import torch
    ctx = {c: v * torch.ones((n,), device=device, dtype=torch.float)
           for c, v in (("logp", round(logp, 2)), ("sascore", round(sa, 2)),
                        ("mol_weight", round(mw, 2)))}
    smiles, _ = sampler.generate(
        context_cols=ctx, context_smi=context_smi, start_smiles=None,
        num_samples=n, max_new_tokens=max_new_tokens, temperature=temperature,
        top_k=top_k, total_gen_steps=int(np.ceil(n / 1000)), return_context=True)
    return smiles


def build_contexts(actives, mode, fragments_path, parents_per_fragment, rng):
    """Return [(context_smi, parent_smiles, kind)] and a diagnostics dict.

    In alert mode every fragment is paired with up to `parents_per_fragment`
    real active molecules that actually contain it, so the property conditions
    come from a real parent and the parents are spread across the splits.
    """
    contexts = []
    diag = {"fragments": 0, "fragments_unparseable": 0, "fragments_zero_match": 0,
            "zero_match_examples": [], "unique_parents": 0}

    active_mols = []
    for s in actives:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            active_mols.append((s, m))

    if mode in ("whole", "mixed"):
        for s, _ in active_mols:
            contexts.append((s, s, "whole"))

    if mode in ("alert", "mixed"):
        if not fragments_path or not os.path.exists(fragments_path):
            raise SystemExit("--fragments is required for alert/mixed mode")
        if os.path.getsize(fragments_path) <= 1:
            raise SystemExit(
                f"{fragments_path} is empty - this assay has no enriched alerts, "
                f"so it cannot have an alert-conditioned arm. That is the expected "
                f"outcome for an alert-depleted assay such as NR-AR; run it with "
                f"--mode whole instead.")
        fr = pd.read_csv(fragments_path)
        if fr.empty:
            raise SystemExit(f"{fragments_path} has a header but no rows.")
        diag["fragments"] = len(fr)
        for _, row in fr.iterrows():
            frag = row["fragment"]
            q = Chem.MolFromSmiles(frag)
            if q is None:
                diag["fragments_unparseable"] += 1
                continue
            matches = [s for s, m in active_mols if m.HasSubstructMatch(q)]
            if not matches:
                # The fragment does not match the molecules it came from, which
                # happens with charge-stripped alert patterns. Using the fragment
                # itself as the parent would make every generation unusable under
                # the clean protocol, so skip it instead.
                diag["fragments_zero_match"] += 1
                if len(diag["zero_match_examples"]) < 8:
                    diag["zero_match_examples"].append(frag)
                continue
            order = rng.permutation(len(matches))
            for j in order[:parents_per_fragment]:
                contexts.append((frag, matches[int(j)], "alert"))

    diag["unique_parents"] = len({p for _, p, _ in contexts})
    return contexts, diag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", default=os.path.join("out", "llama2-M-Full-RSS-Canonical.pt"))
    ap.add_argument("--assay-csv", help="assay data: SMILES,label")
    ap.add_argument("--fragments", help="fragments CSV from step 1 (fragment,alert)")
    ap.add_argument("--mode", choices=["whole", "alert", "mixed"], default="alert")
    ap.add_argument("--out", default="generated.csv")
    ap.add_argument("--per-context", type=int, default=40)
    ap.add_argument("--parents-per-fragment", type=int, default=25,
                    help="alert mode: how many distinct parents to pair with each "
                         "fragment. Raises parent diversity so the clean-protocol "
                         "filter does not drop whole blocks.")
    ap.add_argument("--max-contexts", type=int, default=0,
                    help="0 = no limit")
    ap.add_argument("--target-rows", type=int, default=0,
                    help="how many augmentation rows step 2 will need, used only "
                         "to warn when the projected yield looks too small")
    ap.add_argument("--temperature", type=float, default=0.8)
    ap.add_argument("--top-k", type=int, default=25)
    ap.add_argument("--max-new-tokens", type=int, default=256)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--require-substructure", action="store_true", default=True,
                    help="alert mode: keep only generations that retain the fragment")
    ap.add_argument("--smoke-test", action="store_true",
                    help="generate 5 molecules from one fragment and exit")
    a = ap.parse_args()

    sampler, device = make_sampler(a.ckpt, a.seed)
    print(f"loaded {a.ckpt} on {device}")

    if a.smoke_test:
        frag = "C=CC=O"
        m = Chem.MolFromSmiles(frag)
        lp, sa, mw = props(m)
        out = generate_for(sampler, device, frag, lp, sa, mw, 5,
                           a.temperature, a.top_k, a.max_new_tokens)
        print(f"\ncontext fragment : {frag}  (logp {lp:.2f}, sa {sa:.2f}, mw {mw:.2f})")
        ok = 0
        for s in out:
            mm = Chem.MolFromSmiles(s)
            has = mm is not None and mm.HasSubstructMatch(m)
            ok += bool(has)
            print(f"  valid={mm is not None}  retains_fragment={has}  {s}")
        print(f"\nfragment retained in {ok}/{len(out)}. If this is 0, the model is not "
              f"honouring the SMILES context - stop and check the checkpoint and the "
              f"repo version before running a long job.")
        return

    if not a.assay_csv:
        raise SystemExit("--assay-csv is required unless --smoke-test")

    data = pd.read_csv(a.assay_csv)
    scol = [c for c in data.columns if c.lower() == "smiles"][0]
    minority = int(data["label"].value_counts().idxmin())
    actives = data[data["label"] == minority][scol].tolist()
    known = set()
    for s in data[scol]:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m:
            known.add(Chem.MolToSmiles(m))
    print(f"{a.assay_csv}: minority class {minority}, {len(actives)} actives")

    rng = np.random.default_rng(a.seed)
    contexts, diag = build_contexts(actives, a.mode, a.fragments,
                                    a.parents_per_fragment, rng)
    if a.max_contexts and len(contexts) > a.max_contexts:
        print(f"capping contexts at --max-contexts {a.max_contexts} "
              f"(built {len(contexts)})")
        contexts = contexts[: a.max_contexts]

    if a.mode in ("alert", "mixed"):
        print(f"fragments: {diag['fragments']} total, "
              f"{diag['fragments_unparseable']} unparseable, "
              f"{diag['fragments_zero_match']} matched no active molecule (skipped)")
        if diag["zero_match_examples"]:
            print(f"  skipped examples: {', '.join(diag['zero_match_examples'])}")
    print(f"contexts: {len(contexts)} (mode {a.mode}) over "
          f"{diag['unique_parents']} distinct parents, {a.per_context} molecules each "
          f"-> up to {len(contexts) * a.per_context} generations")
    if not contexts:
        raise SystemExit("no usable contexts were built; nothing to generate")

    if a.target_rows:
        # step 2 keeps only generations whose parent is in the training split,
        # which is about 64 per cent of the data under the default 0.2/0.2 split.
        projected = len(contexts) * a.per_context * 0.64
        print(f"step 2 needs {a.target_rows} rows; projected usable yield before "
              f"validity and duplicate filtering is about {projected:.0f}")
        if projected < 2 * a.target_rows:
            print("  WARNING: that leaves little margin once invalid, duplicate and "
                  "substructure-failing generations are dropped. Raise --per-context "
                  "or --parents-per-fragment.")

    rows, seen = [], set()
    for i, (ctx_smi, parent, kind) in enumerate(contexts, start=1):
        pm = Chem.MolFromSmiles(parent)
        if pm is None:
            continue
        lp, sa, mw = props(pm)
        try:
            gen = generate_for(sampler, device, ctx_smi, lp, sa, mw, a.per_context,
                               a.temperature, a.top_k, a.max_new_tokens)
        except Exception as e:
            print(f"  [{i}/{len(contexts)}] generation failed for {ctx_smi}: {e}")
            continue
        q = Chem.MolFromSmiles(ctx_smi)
        for s in gen:
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            c = Chem.MolToSmiles(m)
            if c in known or c in seen:
                continue
            if kind == "alert" and a.require_substructure and q is not None:
                if not m.HasSubstructMatch(q):
                    continue
            seen.add(c)
            rows.append({"SMILES": c, "label": minority, "parent": parent,
                         "context": ctx_smi, "context_kind": kind,
                         "cond_logp": round(lp, 2), "cond_sascore": round(sa, 2),
                         "cond_mol_weight": round(mw, 2)})
        if i % 10 == 0 or i == len(contexts):
            print(f"  [{i}/{len(contexts)}] kept {len(rows)} unique so far", flush=True)
            pd.DataFrame(rows).to_csv(a.out, index=False)

    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(f"\n{len(rows)} unique generated molecules -> {a.out}")
    if rows:
        df = pd.DataFrame(rows)
        print(f"distinct parents represented: {df['parent'].nunique()}")
        print("columns: SMILES,label,parent,context,context_kind,cond_logp,"
              "cond_sascore,cond_mol_weight")
        print("`parent` is what step 2 uses to keep generations out of the test split.")
        if a.target_rows and len(rows) < a.target_rows:
            print(f"\nWARNING: fewer molecules than the {a.target_rows} rows step 2 "
                  f"needs. Re-run with a larger --per-context before going on.")


if __name__ == "__main__":
    main()
