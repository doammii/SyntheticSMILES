@echo off
REM ===================================================================
REM Stage 3a - generate molecules with LLaMOL.
REM
REM Run this from INSIDE the llamol repository folder, because the script
REM imports the repo's own `sample` module and loads the checkpoint at
REM out\llama2-M-Full-RSS-Canonical.pt
REM
REM TOX21 is set automatically to this repository root.
REM Check a small probe run first (see docs\llamol_probe_result.txt).
REM ===================================================================
setlocal
set TOX21=%~dp0..
set GEN=%TOX21%\generated
if not exist "%GEN%" mkdir "%GEN%"

REM --per-context is sized per assay so that the projected usable yield is
REM about 2.5x the rows stage 2 needs, even on the worst seed. Assays with
REM few alert fragments need a much larger value.
REM --target-rows is only used to print a warning; it changes nothing.

echo ============ SR-ARE (alert) ============
python "%TOX21%\src\03_generate_llamol.py" ^
    --assay-csv "%TOX21%\prepared\assays\SR-ARE.csv" ^
    --fragments "%TOX21%\prepared\fragments\SR-ARE.csv" ^
    --mode alert --per-context 50 --target-rows 1884 ^
    --out "%GEN%\SR-ARE_alert.csv"

echo ============ SR-MMP (alert) ============
python "%TOX21%\src\03_generate_llamol.py" ^
    --assay-csv "%TOX21%\prepared\assays\SR-MMP.csv" ^
    --fragments "%TOX21%\prepared\fragments\SR-MMP.csv" ^
    --mode alert --per-context 80 --target-rows 1836 ^
    --out "%GEN%\SR-MMP_alert.csv"

echo ============ SR-p53 (alert) ============
python "%TOX21%\src\03_generate_llamol.py" ^
    --assay-csv "%TOX21%\prepared\assays\SR-p53.csv" ^
    --fragments "%TOX21%\prepared\fragments\SR-p53.csv" ^
    --mode alert --per-context 180 --target-rows 846 ^
    --out "%GEN%\SR-p53_alert.csv"

echo ============ SR-HSE (alert) ============
python "%TOX21%\src\03_generate_llamol.py" ^
    --assay-csv "%TOX21%\prepared\assays\SR-HSE.csv" ^
    --fragments "%TOX21%\prepared\fragments\SR-HSE.csv" ^
    --mode alert --per-context 200 --target-rows 744 ^
    --out "%GEN%\SR-HSE_alert.csv"

echo ============ NR-PPAR-gamma (alert) ============
python "%TOX21%\src\03_generate_llamol.py" ^
    --assay-csv "%TOX21%\prepared\assays\NR-PPAR-gamma.csv" ^
    --fragments "%TOX21%\prepared\fragments\NR-PPAR-gamma.csv" ^
    --mode alert --per-context 140 --target-rows 372 ^
    --out "%GEN%\NR-PPAR-gamma_alert.csv"

echo ============ NR-AR (whole only - no alert fragments) ============
REM NR-AR is the negative control: BRENK alerts are DEPLETED among its actives
REM (odds ratio 0.70), so prepared\fragments\NR-AR.csv is empty by construction
REM and alert mode is impossible here. Only the whole-molecule arm can run.
python "%TOX21%\src\03_generate_llamol.py" ^
    --assay-csv "%TOX21%\prepared\assays\NR-AR.csv" ^
    --mode whole --per-context 40 --target-rows 616 ^
    --out "%GEN%\NR-AR_whole.csv"

echo.
echo DONE. Send back the console log and the row counts of %GEN%\*.csv
