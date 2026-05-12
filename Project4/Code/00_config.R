# Project 4 simulation configuration
# Run all scripts from the Project4 root directory.

# ==================================================
# Project root directory
# ==================================================
suppressPackageStartupMessages({library(here)})
project_dir <- here() # set working directory location here - use ""
## example here("BIOS6624/Project4")

# ==================================================
# Create folders if they do not exist
# ==================================================

dir.create(here("Code"), showWarnings = FALSE)
dir.create(here("DataProcessed"), showWarnings = FALSE)
dir.create(here("Figures"), showWarnings = FALSE)
dir.create(here("Tables"), showWarnings = FALSE)
dir.create(here("Logs"), showWarnings = FALSE)

set.seed(2522)

n_sim <- 10000
n_values <- c(250, 500)
rho_values <- c(0, 0.35, 0.70)

p <- 20
p_true <- 5
true_vars <- paste0("V", sprintf("%02d", 1:p_true))

beta_true <- c(
  0.5/3,
  1/3,
  1.5/3,
  2/3,
  2.5/3,
  rep(0, 15))
# variable names to match hdrm
names(beta_true) <- paste0("V", sprintf("%02d", 1:p))

alpha_pvalue <- 0.10

alpha_enet <- 0.5
cv_folds <- 10
seed_value <- 2522

