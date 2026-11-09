# SyntheticSMILES

Code and results for **“Revisiting Generative Data Augmentation for Imbalanced Molecular Toxicity Prediction under Leak-Free Evaluation”** (Yong Oh Lee and Do Yeon Kim, 2026).

## Overview

We examine generative SMILES augmentation on six Tox21 toxicity assays (**SR-ARE, SR-MMP, SR-p53, SR-HSE, NR-PPARγ, and NR-AR**) using a public **LLaMol** checkpoint (without fine-tuning) and **ChemBERTa-10M-MLM**. Four conditions—no augmentation, duplication, randomized SMILES, and LLaMol generation—are compared under two protocols:

- **Contaminated:** augmented rows are added before the train/validation/test split.
- **Leak-free:** the original dataset is split first, and augmentation is applied to the training set only.

The study covers **240 training runs** (6 assays × 2 protocols × 4 conditions × 5 seeds), with stratified 64:16:20 train/validation/test splits and both fixed (0.5) and validation-optimized decision thresholds. Augmentation before splitting inflates minority-class AUPRC; under the leak-free protocol, LLaMol augmentation shows no advantage over the simple controls in these six assays. Validation-optimized thresholds reduce many apparent gains observed at a fixed threshold of 0.5.

## Changes from the earlier version

The [earlier preprint](https://doi.org/10.21203/rs.3.rs-6432369/v1) applied augmentation before splitting and incorrectly described LLaMol as fine-tuned. The revised study corrects this description, compares contaminated and leak-free protocols, and adds duplication and randomized-SMILES controls. Its leak-free results supersede the earlier performance conclusions.

## Repository structure

| Path | Contents |
| --- | --- |
| `src/` | Tox21 preprocessing, LLaMol generation, and ChemBERTa experiments. |
| `scripts/` | Windows batch scripts for running the experiments. |
| `analysis/aggregate.py` | Aggregation of raw results across seeds. |
| `data/` | Tox21 input CSV. |
| `prepared/` | Prepared assays, alert-enrichment results, and fragments. |
| `generated/` | Generated molecules used for the LLaMol condition. |
| `results/raw/` | Per-run raw and summary CSVs for all 240 training runs. |
| `docs/` | Supporting experimental notes. |

## Running the experiments

Install Python 3.10, a PyTorch build compatible with your system, and the dependencies in `requirements.txt`:

```bash
pip install -r requirements.txt
```

For generation, clone the [LLaMol repository](https://github.com/Fraunhofer-SCAI/llamol) and place its public `llama2-M-Full-RSS-Canonical.pt` checkpoint in the LLaMol repository's `out/` directory.

1. From the **SyntheticSMILES root** on Windows, prepare the data and run the non-generative conditions:

   ```bat
   python src/01_prepare_tox21.py --outdir prepared --cache data/tox21_original.csv
   call scripts\run_stage2.bat
   ```

2. Optional, to regenerate the molecules in `generated/`. From the **LLaMol repository directory** (substitute the path to your SyntheticSMILES checkout):

   ```bat
   call "C:\path\to\SyntheticSMILES\scripts\run_generate.bat"
   ```

3. Back in the **SyntheticSMILES root**, run the LLaMol condition and aggregate results:

   ```bat
   call scripts\run_classify_llamol.bat
   python analysis/aggregate.py
   ```

The included `results/raw/` files are the reported experimental results, and the CSVs in `generated/` are the generated molecules used to produce them.

## Data and license

- **Tox21:** publicly available via [MoleculeNet](https://moleculenet.org/).
- **Earlier preprint:** https://doi.org/10.21203/rs.3.rs-6432369/v1.
- **License:** repository code is released under the [MIT License](LICENSE). Third-party datasets, models, and dependencies retain their respective terms.
