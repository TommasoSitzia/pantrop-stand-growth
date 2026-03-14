# pantrop-stand-growth

Data and R scripts for analysing age–height and age–DBH growth trajectories in pan-tropical planted forest stands.

**Working manuscript title**  
*A pan-tropical database and synthesis of age–height and age–diameter growth in planted forest stands up to 40 years*

## Overview

This repository contains the datasets and R workflows used to analyse relationships between stand age and two basic tree size attributes in tropical planted forest stands:

- diameter at breast height (DBH)
- total tree height

The repository supports two complementary analyses:

- a **DBH–age analysis**, based on Chapman–Richards growth curves
- a **height–age analysis**, based on a two-component strategy including a baseline Chapman–Richards curve and a separate fast-juvenile power model

The analyses are intended to provide transparent and reproducible default age–size relationships within the age range most strongly supported by the compiled evidence.

## Repository structure

```text
pantrop-stand-growth/
├─ data/
│  ├─ dbhage.xlsx
│  └─ heightage.xlsx
├─ scripts/
│  ├─ dbh_age_analysis.R
│  └─ height_age_analysis.R
├─ outputs/
│  ├─ xlsx/
│  └─ PDF/
│     ├─ DBH/
│     └─ H/
└─ README.md
```

## Input data

The repository expects two Excel files in the `data/` folder:

- `dbhage.xlsx`
- `heightage.xlsx`

In both workbooks, the scripts expect a worksheet named:

- `data`

The scripts are written for **clean plantation-only input datasets**.

### Minimum required columns

For `dbhage.xlsx`:

- `ageexa`
- `ageave`
- `dbhexa`
- `dbhave`
- `taxon`

For `heightage.xlsx`:

- `ageexa`
- `ageave`
- `heightexa`
- `heightave`
- `taxon`

### Additional columns used when available

The scripts can also use additional columns for weighting, grouping, or metadata handling.

For DBH, these may include:

- `dbhsd`, `dbhmin`, `dbhmax`, `notrees`
- `level`
- `country`, `location`
- `idarticle`, `rain`, `temp`, `treeha`, `treespha`, `stand`

For height, these may include:

- `heightsd`, `heightmin`, `heightmax`, `notrees`
- `level`

If some optional columns are absent, the scripts still run as long as the required columns are present.

## Analytical scope

Both analyses use a main modelling window of:

- **age ≤ 40 years**

This is the age range used for the manuscript’s main models and default outputs.

## Script summary

### `scripts/dbh_age_analysis.R`

This script performs the DBH–age analysis.

It:

- reads `data/dbhage.xlsx`
- harmonises age and DBH values by prioritising exact values and otherwise using reported means
- computes optional observation-level weights when uncertainty or sample-size information is available
- restricts the main modelling dataset to `age ≤ 40`
- fits pooled Chapman–Richards models using:
  - unweighted fit
  - capped-weighted fit
  - weighted fit
- uses k-means on taxon mean DBH at ages `≤ 15` to derive DBH growth classes
- fits class-specific Chapman–Richards models for:
  - slow
  - medium
  - fast growth classes
- writes figures, diagnostics, and summary tables

### `scripts/height_age_analysis.R`

This script performs the height–age analysis.

It:

- reads `data/heightage.xlsx`
- harmonises age and height values by prioritising exact values and otherwise using reported means
- computes optional observation-level weights when uncertainty or sample-size information is available
- restricts the main modelling dataset to `age ≤ 40`
- identifies a fast-juvenile subset by clustering taxa on mean height at ages `≤ 15`
- fits a baseline Chapman–Richards model excluding the fast-juvenile group
- fits a fast-juvenile power model over ages `1–10`
- fits a within-window sensitivity Chapman–Richards model including all records
- writes figures, diagnostics, and summary tables

### Note on the height analysis

In the height workflow, the fast-juvenile label is assigned at the **taxon level**, not at the individual-point level. As a result, some observations from fast-juvenile taxa may overlap the baseline cloud. The fast-juvenile power curve itself is fitted only to observations aged `1–10`.

## Outputs

Running the scripts creates local output files in:

- `outputs/PDF/DBH/`
- `outputs/PDF/H/`
- `outputs/xlsx/resultspaper.xlsx`

The two scripts write to the same Excel workbook, but to different worksheet names.

Generated outputs are intended for **local use only** and are **not meant to be tracked in the repository**.

## Required packages

The scripts require the following R packages:

- **dplyr** — data manipulation: filtering, mutating, grouping, counting, and summarising
- **tidyr** — data tidying tools used alongside `dplyr`
- **readxl** — reads the input Excel files
- **ggplot2** — creates the figures and diagnostic plots
- **minpack.lm** — fits non-linear models using the Levenberg–Marquardt algorithm
- **purrr** — applies functions across lists and grouped objects
- **tibble** — creates tidy table outputs
- **openxlsx** — writes results to the Excel workbook
- **gridExtra** — arranges multiple diagnostic plots on a page

Install them with:

```r
install.packages(c(
  "dplyr", "tidyr", "readxl", "ggplot2", "minpack.lm",
  "purrr", "tibble", "openxlsx", "gridExtra"
))
```

## How to run the analyses

Open the repository folder in RStudio and make sure the **working directory is the repository root**.

Then run:

```r
source("scripts/dbh_age_analysis.R")
source("scripts/height_age_analysis.R")
```

The scripts use repository-relative paths via `file.path(...)`, so they must be run from the repository root, not from inside the `scripts/` folder.

## Reproducibility notes

- Both scripts create output folders automatically if needed.
- Both scripts check that the expected input file exists.
- Both scripts check that the expected worksheet name (`data`) is present.
- Random clustering steps use `set.seed(42)`.
- Generated files in `outputs/` are ignored through `.gitignore`.

## Status

This is a working analysis repository associated with a manuscript in preparation. It is intended primarily for reproducibility, internal verification, and future public release alongside the paper.

## Contact

Tommaso Sitzia  
Department of Land, Environment, Agriculture and Forestry  
University of Padova
tommaso.sitzia@unipd.it
