@echo off
REM Run from the repository root with the conda env active (see README).
REM Order: 1) prepare  2) baselines  3) generate (needs llamol repo)  4) llamol condition
setlocal
python src\01_prepare_tox21.py --outdir prepared --cache data\tox21_original.csv
call scripts\run_stage2.bat
REM Step 3 must be run inside a clone of the llamol repository, see scripts\run_generate.bat
call scripts\run_classify_llamol.bat
python analysis\aggregate.py
