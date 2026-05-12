#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(here)
})

source(here("Code", "00_config.R"))

options(scipen = 999)

method_summary <- readRDS(file.path(
  project_dir, "DataProcessed", "method_summary.rds"))
variable_summary <- readRDS(file.path(
  project_dir, "DataProcessed", "variable_summary.rds"))
mcse_inference_table <- readRDS(file.path(
  project_dir, "DataProcessed", "mcse_inference_table.rds"))

method_labels <- c(
  backward_pvalue = "Backward (p-value)",
  backward_AIC = "AIC",
  backward_BIC = "BIC",
  lasso_lambda_min = "LASSO (min)",
  lasso_lambda_1se = "LASSO (1se)",
  elastic_net_lambda_min = "ElNet (min)",
  elastic_net_lambda_1se = "ElNet (1se)")

# Figure 1: Bias heatmap
bias_heatmap_data <- variable_summary %>%
  filter(beta_true != 0) %>%
  mutate(Method = recode(method, !!!method_labels))

fig_bias <- ggplot(bias_heatmap_data, 
                   aes(x = variable, y = Method, fill = mean_bias)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_bias, 2)), size = 5) +
  facet_grid(n ~ rho) +
  scale_fill_gradient2(
    breaks = c(-0.09, 0, 0.02),
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "Mean Bias") +
  labs(x = "True Predictor", y = "Method", 
       title = "Mean Bias for True Predictors") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom")

ggsave(here("Figures", "figure1_bias_heatmap.png"), 
       fig_bias, width = 12, height = 8, dpi = 300)

# Figure 3: Coverage heatmap
coverage_data <- variable_summary %>%
  filter(beta_true != 0) %>%
  mutate(Method = recode(method, !!!method_labels))

fig_coverage <- ggplot(coverage_data, 
                       aes(x = variable, y = Method, fill = coverage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(coverage, 3), color = coverage < 0.8), size = 5) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  facet_grid(n ~ rho) +
  scale_fill_gradientn(
    colours = c("#67001F", "#B2182B", "#F4A582", "#FEEBE2", "white", 
                "#D9F0D3", "#1B7837"),
    values = scales::rescale(c(0.20, 0.60, 0.85, 0.93, 0.95, 0.97, 1.00), 
                             from = c(0.20, 1.00)),
    limits = c(0.20, 1.00),
    name = "Coverage") +
  labs(x = "True Predictor", y = "Method", 
       title = "Coverage of 95% CIs for True Predictors") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom")

ggsave(here("Figures", "figure3_coverage_heatmap.png"), 
       fig_coverage, width = 12, height = 8, dpi = 300)

# Supplemental Figure 1: Selection probability heatmap
selection_heatmap_data <- variable_summary %>%
  filter(beta_true != 0) %>%
  mutate(Method = recode(method, !!!method_labels))

fig_selection <- 
  ggplot(selection_heatmap_data, 
         aes(x = variable, y = Method, fill = selection_probability)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(selection_probability, 3), 
                color = selection_probability < 0.75), size = 5) +
  facet_grid(n ~ rho) +
  scale_fill_gradientn(
    colours = c("#B2182B", "#EF8A62", "#FEE0D2", "white", "#D9F0D3", "#1B7837"),
    values = scales::rescale(c(0.20, 0.50, 0.75, 0.90, 0.97, 1.00), 
                             from = c(0.20, 1.00)),
    limits = c(0.20, 1.00),
    name = "Selection\nProbability") +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  labs(x = "True Predictor", y = "Method", 
       title = "Selection Probability for True Predictors") +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom")


ggsave(here("Figures", "supplement_selection_probability.png"), 
       fig_selection, width = 12, height = 8, dpi = 300)

# Supplemental Figure 2: MCSE selection performance
mcse_heatmap_data <- method_summary %>%
  mutate(Method = recode(method, !!!method_labels)) %>%
  select(n, rho, Method, `TPR` = mcse_sensitivity, 
         `FPR` = mcse_false_positive_rate, `Model Size` = mcse_model_size) %>%
  pivot_longer(cols = c(`TPR`, `FPR`, `Model Size`), names_to = "metric", 
               values_to = "mcse")

fig_mcse_selection <- ggplot(mcse_heatmap_data, 
                             aes(x = factor(rho), y = Method, fill = mcse)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mcse, 4), color = mcse > 0.02), size = 5) +
  facet_grid(metric ~ n) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "MCSE") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), 
                     guide = "none") +
  labs(x = "Correlation (rho)", y = "Method", 
       title = "Monte Carlo Standard Errors for Selection Performance") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("Figures", "supplement_mcse_selection.png"), 
       fig_mcse_selection, width = 12, height = 8, dpi = 300)

# Supplemental Figure 3: MCSE inferential performance
mcse_inference_heatmap <- mcse_inference_table %>%
  mutate(Method = recode(method, !!!method_labels)) %>%
  pivot_longer(cols = c(mcse_bias, mcse_coverage, mcse_type2, mcse_type1), 
               names_to = "metric", values_to = "mcse") %>%
  mutate(metric = recode(metric, mcse_bias = "Bias", 
                         mcse_coverage = "Coverage", 
                         mcse_type2 = "Type II", 
                         mcse_type1 = "Type I"))

fig_mcse_inference <- ggplot(mcse_inference_heatmap, 
                             aes(x = factor(rho), y = Method, fill = mcse)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mcse, 4), color = mcse > 0.002), size = 5) +
  facet_grid(metric ~ n) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "MCSE") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"),
                     guide = "none") +
  labs(x = "Correlation (rho)", y = "Method", 
       title = "Monte Carlo Standard Errors for Inferential Performance") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("Figures", "supplement_mcse_inference.png"), 
       fig_mcse_inference, width = 12, height = 8, dpi = 300)

message("Saved figures to: ", file.path(project_dir, "Figures"))
