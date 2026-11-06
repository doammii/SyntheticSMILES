@echo off
REM Stage 2 - classifier experiments. Activate the env first: conda activate tox21
REM Six assays spanning alert enrichment and imbalance, both protocols.
setlocal
for %%A in (SR-ARE SR-HSE SR-p53 NR-PPAR-gamma SR-MMP NR-AR) do (
  for %%P in (clean contaminated) do (
    echo ============ %%A : %%P ============
    python src/02_run_experiments.py --data prepared/assays/%%A.csv ^
        --outdir results/raw/%%A_%%P --protocol %%P ^
        --conditions baseline dup randsmi --seeds 0 1 2 3 4
  )
)
echo.
echo DONE. Results are in results\raw\*
