#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(gt)
  library(here)
})

source(here("Code", "00_config.R"))

method_summary <- readRDS(
  here("DataProcessed", "method_summary.rds"))

variable_summary <- readRDS(
  here("DataProcessed", "variable_summary.rds"))

# This script saves gt tables as HTML.

method_labels <- c(
  backward_pvalue = "Backward (p-value)",
  backward_AIC = "AIC",
  backward_BIC = "BIC",
  lasso_lambda_min = "LASSO (λₘᵢₙ)",
  lasso_lambda_1se = "LASSO (λ₁ₛₑ)",
  elastic_net_lambda_min = "ElNet (λₘᵢₙ)",
  elastic_net_lambda_1se = "ElNet (λ₁ₛₑ)")

rho0_cols <- c("rho0_n250_TPR", "rho0_n250_FPR", "rho0_n250_Size", 
               "rho0_n500_TPR", "rho0_n500_FPR", "rho0_n500_Size")
rho035_cols <- c("rho035_n250_TPR", "rho035_n250_FPR", "rho035_n250_Size", 
                 "rho035_n500_TPR", "rho035_n500_FPR", "rho035_n500_Size")
rho07_cols <- c("rho07_n250_TPR", "rho07_n250_FPR", "rho07_n250_Size", 
                "rho07_n500_TPR", "rho07_n500_FPR", "rho07_n500_Size")

selection_table_wide <- method_summary %>%
  mutate(
    Method = recode(method, !!!method_labels),
    rho_label = case_when(
      rho == 0 ~ "rho0", 
      rho == 0.35 ~ "rho035", 
      rho == 0.7 ~ "rho07")) %>%
  select(Method, rho_label, n, 
         TPR = mean_sensitivity, 
         FPR = mean_false_positive_rate, 
         Size = mean_model_size) %>%
  pivot_wider(names_from = c(rho_label, n), 
              values_from = c(TPR, FPR, Size), 
              names_glue = "{rho_label}_n{n}_{.value}") %>%
  select(Method, all_of(rho0_cols), all_of(rho035_cols), all_of(rho07_cols))

selection_gt <- selection_table_wide %>%
  gt(rowname_col = "Method") %>%
  tab_header(title = "Variable Selection Performance") %>%
  tab_spanner(label = "n = 250", id = "rho0_n250", 
              columns = c(rho0_n250_TPR, rho0_n250_FPR, rho0_n250_Size)) %>%
  tab_spanner(label = "n = 500", id = "rho0_n500", 
              columns = c(rho0_n500_TPR, rho0_n500_FPR, rho0_n500_Size)) %>%
  tab_spanner(label = "n = 250", id = "rho035_n250", 
              columns = c(rho035_n250_TPR, rho035_n250_FPR, rho035_n250_Size)) %>%
  tab_spanner(label = "n = 500", id = "rho035_n500", 
              columns = c(rho035_n500_TPR, rho035_n500_FPR, rho035_n500_Size)) %>%
  tab_spanner(label = "n = 250", id = "rho07_n250", 
              columns = c(rho07_n250_TPR, rho07_n250_FPR, rho07_n250_Size)) %>%
  tab_spanner(label = "n = 500", id = "rho07_n500", 
              columns = c(rho07_n500_TPR, rho07_n500_FPR, rho07_n500_Size)) %>%
  tab_spanner(label = "ρ = 0", id = "rho0", 
              columns = all_of(rho0_cols)) %>%
  tab_spanner(label = "ρ = 0.35", id = "rho035", 
              columns = all_of(rho035_cols)) %>%
  tab_spanner(label = "ρ = 0.7", id = "rho07", 
              columns = all_of(rho07_cols)) %>%
  cols_label(.list = setNames(rep(c("TPR", "FPR", "Size"), 6), 
                              names(selection_table_wide)[-1])) %>%
  fmt_number(columns = contains("TPR") | contains("FPR"), decimals = 2) %>%
  fmt_number(columns = contains("Size"), decimals = 1)

gtsave(selection_gt, here("Tables", "table1_selection_performance.html"))

error_table <- variable_summary %>%
  group_by(n, rho, method) %>%
  summarise(
    type1 = mean(type1_error_rate[beta_true == 0], na.rm = TRUE),
    type2 = mean(type2_error_rate[beta_true != 0], na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    Method = recode(method, !!!method_labels),
    rho_label = case_when(
      rho == 0 ~ "rho0", 
      rho == 0.35 ~ "rho035", 
      rho == 0.7 ~ "rho07")) %>%
  select(Method, rho_label, n, type1, type2) %>%
  pivot_wider(names_from = c(rho_label, n), 
              values_from = c(type1, type2), 
              names_glue = "{rho_label}_n{n}_{.value}") %>%
  select(
    Method,
    rho0_n250_type1, rho0_n250_type2, rho0_n500_type1, rho0_n500_type2,
    rho035_n250_type1, rho035_n250_type2, rho035_n500_type1, rho035_n500_type2,
    rho07_n250_type1, rho07_n250_type2, rho07_n500_type1, rho07_n500_type2)

error_gt <- error_table %>%
  gt(rowname_col = "Method") %>%
  tab_header(title = "Type I and Type II Error Rates") %>%
  tab_spanner(label = "n = 250", id = "rho0_n250", 
              columns = c(rho0_n250_type1, rho0_n250_type2)) %>%
  tab_spanner(label = "n = 500", id = "rho0_n500", 
              columns = c(rho0_n500_type1, rho0_n500_type2)) %>%
  tab_spanner(label = "n = 250", id = "rho035_n250", 
              columns = c(rho035_n250_type1, rho035_n250_type2)) %>%
  tab_spanner(label = "n = 500", id = "rho035_n500", 
              columns = c(rho035_n500_type1, rho035_n500_type2)) %>%
  tab_spanner(label = "n = 250", id = "rho07_n250", 
              columns = c(rho07_n250_type1, rho07_n250_type2)) %>%
  tab_spanner(label = "n = 500", id = "rho07_n500", 
              columns = c(rho07_n500_type1, rho07_n500_type2)) %>%
  tab_spanner(label = "ρ = 0", id = "rho0", 
              columns = c(rho0_n250_type1, rho0_n250_type2, rho0_n500_type1, 
                          rho0_n500_type2)) %>%
  tab_spanner(label = "ρ = 0.35", id = "rho035", 
              columns = c(rho035_n250_type1, rho035_n250_type2, 
                          rho035_n500_type1, rho035_n500_type2)) %>%
  tab_spanner(label = "ρ = 0.7", id = "rho07", 
              columns = c(rho07_n250_type1, rho07_n250_type2, 
                          rho07_n500_type1, rho07_n500_type2)) %>%
  cols_label(
    rho0_n250_type1 = "Type I", rho0_n250_type2 = "Type II", 
    rho0_n500_type1 = "Type I", rho0_n500_type2 = "Type II",
    rho035_n250_type1 = "Type I", rho035_n250_type2 = "Type II", 
    rho035_n500_type1 = "Type I", rho035_n500_type2 = "Type II",
    rho07_n250_type1 = "Type I", rho07_n250_type2 = "Type II", 
    rho07_n500_type1 = "Type I", rho07_n500_type2 = "Type II") %>%
  fmt_number(columns = starts_with("rho"), decimals = 3)

gtsave(error_gt, here("Tables", "table2_type1_type2_errors.html"))

message("Saved tables to: ", file.path(project_dir, "Tables"))
