# Project 4 Reproducibility Instructions

All analyses for this project were conducted using scripted workflows in R. 
The simulation, summary results, tables, and figures can be reproduced directly 
from the provided code.

---

# Software

The project was built using:

- R (version 4.5.3 "Reassured Reassurer")

## Required R Packages

The following packages are required for the scripts and/or R Markdown document:

```r
library(hdrm)
library(glmnet)
library(MASS)
library(dplyr)
library(tibble)
library(future)
library(furrr)
library(progress)
library(progressr)
library(ggplot2)
library(future.apply)
library(gt)
library(gtsummary)
library(tidyr)
library(here)
```

---

# Project Structure

```text
Project4/
├── Code/
│   ├── 00_config.R
│   ├── 01_helpers.R
│   ├── 02_run_simulation.R
│   ├── 03_summarize_results.R
│   ├── 04_make_figures.R
│   └── 05_make_tables.R
│
├── DataProcessed/
├── Figures/
├── Tables/
│
└── run_all.R
```

Directories such as `DataProcessed/`, `Figures/`, and `Tables/` are created 
automatically from the scripts if they do not already exist.

---

# Running the Entire Project

The entire simulation study can be reproduced by running:

```r
source("run_all.R")
```

Alternatively, from a terminal:

```bash
Rscript run_all.R
```

This will sequentially:

1. Run all simulation scenarios
2. Save simulation outputs
3. Generate summary statistics
4. Create figures
5. Create tables

---

# Running Individual Components

Scripts may also be run individually:

```r
library(here)

source(here("Code", "00_config.R"))
source(here("Code", "01_helpers.R"))
source(here("Code", "02_run_simulation.R"))
source(here("Code", "03_summarize_results.R"))
source(here("Code", "04_make_figures.R"))
source(here("Code", "05_make_tables.R"))
```

---

# Reproducing the Report

The report can also be reproduced by rendering the R Markdown document:

```r
rmarkdown::render("Code/Project4_Sim.Rmd")
```

or from a terminal:

```bash
Rscript -e "rmarkdown::render('Code/Project4_Sim.Rmd')"
```

or by simply running the .RMD within Rstudio

---

# Notes

- Be sure to check the directory settings in `00_config.R` if running:

```bash
Rscript run_all.R
```

- Alpine-compatible scripts are provided. A `p4_sim.batch.sh` file is included 
for running the simulation on Alpine.

- If using Alpine, review the worker/core settings in:

```text
Code/02_run_simulation.R
```

(lines 101–104).

- Simulation outputs are saved as `.rds` files within the `DataProcessed/` 
directory when running the standalone `.R` scripts instead of the R Markdown 
document.

- Any additional questions, feel free to reach out.
