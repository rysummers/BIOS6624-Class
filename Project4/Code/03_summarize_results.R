#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(here)
})

source(here("Code", "00_config.R"))
source(here("Code", "01_helpers.R"))

variable_results <- readRDS(file.path(
  project_dir, "DataProcessed", "variable_results.rds"))
scenario_results <- readRDS(file.path(
  project_dir, "DataProcessed", "scenario_results.rds"))

method_summary <- scenario_results %>%
  group_by(n, rho, method) %>%
  summarise(
    mean_TP = mean(TP, na.rm = TRUE),
    mcse_TP = mcse_mean(TP),
    mean_FP = mean(FP, na.rm = TRUE),
    mcse_FP = mcse_mean(FP),
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    mcse_sensitivity = mcse_mean(sensitivity),
    mean_specificity = mean(specificity, na.rm = TRUE),
    mcse_specificity = mcse_mean(specificity),
    mean_false_positive_rate = mean(false_positive_rate, na.rm = TRUE),
    mcse_false_positive_rate = mcse_mean(false_positive_rate),
    mean_model_size = mean(model_size, na.rm = TRUE),
    mcse_model_size = mcse_mean(model_size),
    .groups = "drop")

variable_summary <- variable_results %>%
  group_by(n, rho, method, variable, beta_true) %>%
  summarise(
    selection_probability = mean(selected),
    mcse_selection_probability = mcse_prop(selected),
    mean_bias = mean(bias, na.rm = TRUE),
    mcse_bias = mcse_mean(bias),
    mse = mean(bias^2, na.rm = TRUE),
    mcse_mse = mcse_mean(bias^2),
    coverage = mean(covered, na.rm = TRUE),
    mcse_coverage = mcse_prop(covered),
    type1_error_rate = ifelse(unique(beta_true) == 0, 
                              mean(type1_error), NA_real_),
    mcse_type1_error = ifelse(unique(beta_true) == 0, 
                              mcse_prop(type1_error), NA_real_),
    type2_error_rate = ifelse(unique(beta_true) != 0, 
                              mean(type2_error), NA_real_),
    mcse_type2_error = ifelse(unique(beta_true) != 0, 
                              mcse_prop(type2_error), NA_real_),
    .groups = "drop")

mcse_inference_table <- variable_summary %>%
  group_by(n, rho, method) %>%
  summarise(
    mcse_bias = mean(mcse_bias[beta_true != 0], na.rm = TRUE),
    mcse_coverage = mean(mcse_coverage[beta_true != 0], na.rm = TRUE),
    mcse_type1 = mean(mcse_type1_error[beta_true == 0], na.rm = TRUE),
    mcse_type2 = mean(mcse_type2_error[beta_true != 0], na.rm = TRUE),
    .groups = "drop")

saveRDS(method_summary, file.path(
  project_dir, "DataProcessed", "method_summary.rds"))
saveRDS(variable_summary, file.path(
  project_dir, "DataProcessed", "variable_summary.rds"))
saveRDS(mcse_inference_table, file.path(
  project_dir, "DataProcessed", "mcse_inference_table.rds"))

message("Saved summaries to: ", file.path(project_dir, "DataProcessed"))
