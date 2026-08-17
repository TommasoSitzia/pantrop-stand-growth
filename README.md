# pantrop-stand-growth

Data and R scripts for analysing age–height and age–DBH growth trajectories in pan-tropical planted forest stands.

## Working manuscript title

**A pan-tropical database and synthesis of age–height and age–diameter growth in planted forest stands up to 40 years**

## Overview

This repository contains the datasets and R workflows used to analyse relationships between stand age and two basic tree-size attributes in tropical planted forest stands:

- diameter at breast height (DBH)
- total tree height

The repository supports two complementary analyses:

- a **DBH–age analysis** comprising a pooled Chapman–Richards reference curve, Chapman–Richards curves for empirically defined slow and medium growth classes, and a non-asymptotic power model for the fast DBH class;
- a **height–age analysis** based on a two-component strategy comprising a baseline Chapman–Richards curve and a separate fast-juvenile power model.

The analyses are intended to provide transparent and reproducible general empirical reference relationships and alternative growth scenarios within the age range most strongly supported by the compiled evidence.

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
│  │  ├─ dbh_results.xlsx
│  │  └─ height_results.xlsx
│  └─ PDF/
│     ├─ DBH/
│     └─ H/
└─ README.md
````

## Input data

The repository expects two Excel files in the `data/` folder:

* `dbhage.xlsx`
* `heightage.xlsx`

In both workbooks, the scripts expect a worksheet named:

```text
data
```

The scripts are written for clean plantation-only input datasets.

### Minimum required columns

For `dbhage.xlsx`:

```text
ageexa
ageave
dbhexa
dbhave
taxon
```

For `heightage.xlsx`:

```text
ageexa
ageave
heightexa
heightave
taxon
```

### Additional columns used when available

The scripts use additional columns for weighting, study-aware analyses, reporting-level sensitivity, and metadata handling.

For DBH, these may include:

```text
idarticle
dbhsd
dbhmin
dbhmax
notrees
level
country
location
rain
temp
treeha
treespha
stand
```

For height, these may include:

```text
idarticleh
heightsd
heightmin
heightmax
notrees
level
country
location
```

The study identifiers `idarticle` and `idarticleh` are required when the study-balanced and leave-one-study-out analyses are enabled.

## Analytical scope

Both analyses use a main modelling window of:

```text
age ≤ 40 years
```

This is the age range used for the manuscript's main models and associated sensitivity analyses.

The observations are concentrated substantially earlier within this interval, so the fitted relationships should be interpreted primarily as empirical references for young and mid-rotation planted stands rather than as long-term extrapolation models.

## Script summary

### `scripts/dbh_age_analysis.R`

This script performs the DBH–age analysis.

It:

* reads `data/dbhage.xlsx`;
* harmonises age and DBH values by prioritising exact values and otherwise using reported means;
* computes optional observation-level weights when uncertainty or sample-size information is available;
* restricts the main modelling dataset to age ≤ 40 years;
* fits pooled Chapman–Richards models using:

  * an unweighted fit used as the principal general DBH reference;
  * a capped-weighted fit;
  * a fully weighted fit;
  * a study-balanced fit;
* uses k-means clustering on taxon mean DBH at ages ≤ 15 years to derive empirical slow, medium, and fast growth classes;
* fits Chapman–Richards relationships to the slow and medium growth classes;
* compares Chapman–Richards, power, and log-linear alternatives for the fast DBH class;
* uses a non-asymptotic power relationship as the selected fast-class representation within its observed age range;
* evaluates parameter uncertainty;
* evaluates the distribution and influence of observation-level weights;
* quantifies study dominance;
* performs leave-one-study-out validation;
* evaluates sensitivity of the growth-class classification to alternative thresholds;
* evaluates sensitivity to reporting level;
* writes figures, diagnostic plots, predictions, and summary tables.

The fast DBH power relationship is intended for interpretation within the observed age domain and is not used for unsupported long-term extrapolation.

### `scripts/height_age_analysis.R`

This script performs the height–age analysis.

It:

* reads `data/heightage.xlsx`;
* harmonises age and height values by prioritising exact values and otherwise using reported means;
* computes optional observation-level weights when uncertainty or sample-size information is available;
* restricts the main modelling dataset to age ≤ 40 years;
* identifies a fast-juvenile component by clustering taxa according to mean height at ages ≤ 15 years;
* assigns the fast-juvenile label at the taxon level;
* fits a baseline Chapman–Richards model after excluding records belonging to fast-juvenile taxa;
* fits a separate fast-juvenile power model using observations aged 1–10 years;
* fits a Chapman–Richards model including all height records as a within-window sensitivity analysis;
* compares unweighted, capped-weighted, weighted, and study-balanced baseline fits;
* evaluates parameter uncertainty;
* quantifies study dominance;
* performs leave-one-study-out validation for the baseline and fast-juvenile relationships;
* evaluates classification sensitivity under alternative age and minimum-observation thresholds;
* evaluates sensitivity of the fast-juvenile power model to alternative fitting windows;
* compares the operational baseline with a baseline restricted to explicitly classified taxa;
* evaluates sensitivity to reporting level;
* writes figures, diagnostic plots, predictions, and summary tables.

## Note on the height analysis

The fast-juvenile label is assigned at the **taxon level**, not independently to individual observations. Consequently, some observations from taxa classified as fast-juvenile may overlap the baseline cloud.

The fast-juvenile power relationship itself is fitted only to observations aged 1–10 years and is interpreted as an empirical upper juvenile-growth scenario rather than as a mature-height trajectory.

Taxa that do not meet the minimum requirements for clustering are retained in the operational baseline unless they belong to a taxon classified as fast-juvenile.

## Outputs

Running the scripts creates local output files under:

```text
outputs/PDF/DBH/
outputs/PDF/H/
outputs/xlsx/dbh_results.xlsx
outputs/xlsx/height_results.xlsx
```

The PDF folders contain manuscript figures, sensitivity figures, and model-diagnostic plots.

The Excel workbooks contain model coefficients, parameter uncertainty, classification summaries, study-dominance results, leave-one-study-out validation results, sensitivity analyses, and associated predictions.

Generated outputs are intended for local analysis and verification and do not need to be tracked in the repository.

## Required packages

The scripts require the following R packages:

* `dplyr` — data manipulation, filtering, grouping, and summarisation
* `tidyr` — data tidying
* `readxl` — reading input Excel files
* `ggplot2` — figures and diagnostic plots
* `minpack.lm` — nonlinear least-squares fitting using the Levenberg–Marquardt algorithm
* `purrr` — functional iteration across models and scenarios
* `tibble` — tidy table outputs
* `openxlsx` — writing Excel output workbooks
* `gridExtra` — arranging diagnostic plots

Install them with:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "readxl",
  "ggplot2",
  "minpack.lm",
  "purrr",
  "tibble",
  "openxlsx",
  "gridExtra"
))
```

## How to run the analyses

Clone or download the repository and set the repository root as the current working directory in R or RStudio.

Then run:

```r
source("scripts/dbh_age_analysis.R")
source("scripts/height_age_analysis.R")
```

The scripts use repository-relative paths and therefore should be run from the repository root.

Expected input paths are:

```text
data/dbhage.xlsx
data/heightage.xlsx
```

Outputs are written automatically to the `outputs/` directory.

## Reproducibility notes

* Both scripts create required output folders automatically if they do not already exist.
* Both scripts check that the expected input file exists.
* Both scripts check that the expected worksheet named `data` is present.
* Required input columns are checked before model fitting.
* Random clustering procedures use `set.seed(42)`.
* Principal models and sensitivity analyses use the same harmonised age ≤ 40-year modelling datasets.
* Study-balanced fits reduce the influence of studies contributing disproportionately large numbers of observations.
* Leave-one-study-out validation removes entire source studies at each iteration to evaluate between-study transferability without record-level leakage.
* Classification-sensitivity analyses evaluate whether empirical growth groups are stable under alternative eligibility thresholds and early-age windows.
* Generated files in `outputs/` can be excluded from version control through `.gitignore`.

## Interpretation

The fitted functions are intended as broad empirical reference trajectories for situations in which stand age is known but a suitable local species- or site-specific growth equation is unavailable.

They are not intended to replace locally calibrated growth-and-yield models where adequate local data exist.

The pooled DBH and baseline height relationships provide general references, while the slow and medium DBH classes and the fast DBH and fast-juvenile height relationships provide alternative empirical growth scenarios within their supported domains.

## Status

This repository accompanies a manuscript under revision.

It is intended for:

* reproducibility of the published analyses;
* transparent inspection of the underlying data and analytical workflow;
* independent verification;
* future extension as additional eligible observations become available.

The database and scripts may be updated before the final archived release associated with the published article.

## Contact

**Tommaso Sitzia**
Department of Land, Environment, Agriculture and Forestry
University of Padova
[tommaso.sitzia@unipd.it](mailto:tommaso.sitzia@unipd.it)

