@echo off
REM ===================================================================
REM Stage 3b - add the llamol condition to the stage 2 experiments.
REM
REM Run AFTER generation, from the tox21_augmentation folder, with the
REM tox21 environment active:  conda activate tox21
REM
REM Why --n-add is written out explicitly for every assay:
REM   Stage 2 defaults the number of added rows to 2 x the minority count
REM   when no synthetic file is given - that is what produced the runs of
REM   2026-10-06. As soon as --synthetic is passed, the default becomes
REM   len(synthetic) instead, which would change the augmentation volume
REM   and make the new numbers incomparable with those runs. Pinning
REM   --n-add to the values below keeps every condition volume-matched.
REM
REM Results go to NEW folders (out\<assay>_<protocol>_llamol) so the
REM 2026-10-06 results are left untouched.
REM ===================================================================
setlocal

for %%P in (clean contaminated) do (
  echo ============ SR-ARE : %%P : llamol ============
  python src/02_run_experiments.py --data prepared/assays/SR-ARE.csv ^
      --synthetic generated/SR-ARE_alert.csv --n-add 1884 ^
      --outdir results/raw/SR-ARE_%%P_llamol --protocol %%P ^
      --conditions llamol --seeds 0 1 2 3 4

  echo ============ SR-MMP : %%P : llamol ============
  python src/02_run_experiments.py --data prepared/assays/SR-MMP.csv ^
      --synthetic generated/SR-MMP_alert.csv --n-add 1836 ^
      --outdir results/raw/SR-MMP_%%P_llamol --protocol %%P ^
      --conditions llamol --seeds 0 1 2 3 4

  echo ============ SR-p53 : %%P : llamol ============
  python src/02_run_experiments.py --data prepared/assays/SR-p53.csv ^
      --synthetic generated/SR-p53_alert.csv --n-add 846 ^
      --outdir results/raw/SR-p53_%%P_llamol --protocol %%P ^
      --conditions llamol --seeds 0 1 2 3 4

  echo ============ SR-HSE : %%P : llamol ============
  python src/02_run_experiments.py --data prepared/assays/SR-HSE.csv ^
      --synthetic generated/SR-HSE_alert.csv --n-add 744 ^
      --outdir results/raw/SR-HSE_%%P_llamol --protocol %%P ^
      --conditions llamol --seeds 0 1 2 3 4

  echo ============ NR-PPAR-gamma : %%P : llamol ============
  python src/02_run_experiments.py --data prepared/assays/NR-PPAR-gamma.csv ^
      --synthetic generated/NR-PPAR-gamma_alert.csv --n-add 372 ^
      --outdir results/raw/NR-PPAR-gamma_%%P_llamol --protocol %%P ^
      --conditions llamol --seeds 0 1 2 3 4

  echo ============ NR-AR : %%P : llamol, whole-molecule context ============
  python src/02_run_experiments.py --data prepared/assays/NR-AR.csv ^
      --synthetic generated/NR-AR_whole.csv --n-add 616 ^
      --outdir results/raw/NR-AR_%%P_llamol --protocol %%P ^
      --conditions llamol --seeds 0 1 2 3 4
)

echo.
echo DONE. Send back out\*_llamol\results_summary.csv, out\*_llamol\results_raw.csv
echo and the console log. Watch for any line starting with "WARNING:" - it means
echo a seed did not have enough generated molecules to match the other conditions.
