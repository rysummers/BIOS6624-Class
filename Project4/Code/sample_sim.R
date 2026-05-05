
# Project 4 Simulation Skeleton ============================================
# Variable Selection Algorithms in Linear Regression

# Install if needed:
# install.packages(c("glmnet", "MASS", "dplyr", "purrr", "tibble"))
# remotes::install_github("pbreheny/hdrm")

library(hdrm)
library(glmnet)
library(MASS)
library(dplyr)
library(tibble)
library(future)
library(furrr)
library(progress)
library(progressr)
library(qs)
library(ggplot2)
library(future.apply)
library(gt)
library(gtsummary)

# load objects
scenario_results <- readRDS(file.path(
  "/Users/ryan_summers/Library/Mobile\ Documents/com~apple~CloudDocs",
  "GitHub/BIOS6624-Class/Project4/DataProcessed/scenario_results.rds"))
variable_results <- readRDS(file.path(
  "/Users/ryan_summers/Library/Mobile\ Documents/com~apple~CloudDocs",
  "GitHub/BIOS6624-Class/Project4/DataProcessed/variable_results.rds"))

scenario_results2 <- readRDS("/Users/ryan_summers/Desktop/scenario_results.rds")
variable_results2 <- readRDS("/Users/ryan_summers/Desktop/variable_results.rds")

## Simulation settings ==================================================== 

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

names(beta_true) <- paste0("V", sprintf("%02d", 1:p))

alpha_pvalue <- 0.10 # change this to 0.10 - stated to be similar to AIC @ 0.157
alpha_enet <- 0.5 # perhaps test multiple alphas??? c(0.25, 0.5, 0.75)
cv_folds <- 10


## Helper Functions ==================================================== 

backward_pvalue <- function(full_model, alpha_remove = 0.157) {
  current_model <- full_model
  
  repeat {
    pvals <- summary(current_model)$coefficients[-1, 4]
    
    if (length(pvals) == 0) break
    
    max_p <- max(pvals, na.rm = TRUE)
    
    if (max_p <= alpha_remove) break
    
    remove_var <- names(which.max(pvals))
    current_terms <- attr(terms(current_model), "term.labels")
    updated_terms <- setdiff(current_terms, remove_var)
    
    if (length(updated_terms) == 0) {
      current_model <- lm(
        formula = as.formula("y ~ 1"),
        data = model.frame(current_model))
      break
      }
    
    new_formula <- as.formula(
      paste("y ~", paste(updated_terms, collapse = " + ")))
    
    current_model <- lm(new_formula, data = model.frame(current_model))
    }
  
  current_model
  }

get_selected_lm <- function(model) {
  setdiff(names(coef(model)), "(Intercept)")
  }

get_selected_glmnet <- function(cvfit, s_value) {
  b <- coef(cvfit, s = s_value)
  selected <- rownames(b)[as.vector(b != 0)]
  setdiff(selected, "(Intercept)")
  }

post_selection_refit <- function(y, X, selected_vars) {
  df <- data.frame(y = y, X)
  
  if (length(selected_vars) == 0) {
    return(tibble(
      variable = character(),
      estimate = numeric(),
      p_value = numeric(),
      ci_low = numeric(),
      ci_high = numeric()))
    }
  
  form <- as.formula(paste("y ~", paste(selected_vars, collapse = " + ")))
  fit <- lm(form, data = df)
  
  coefs <- summary(fit)$coefficients
  cis <- confint(fit)
  
  tibble(
    variable = rownames(coefs)[-1],
    estimate = coefs[-1, "Estimate"],
    p_value = coefs[-1, "Pr(>|t|)"],
    ci_low = cis[-1, 1],
    ci_high = cis[-1, 2])
  }

evaluate_method <- function(
    method_name, 
    selected_vars, 
    y, 
    X, 
    beta_true, 
    alpha_test = 0.05) {
  
  all_vars <- names(beta_true)
  true_nonzero <- names(beta_true[beta_true != 0])
  true_null <- names(beta_true[beta_true == 0])
  
  refit <- post_selection_refit(y, X, selected_vars)
  
  variable_level <- tibble(variable = all_vars) %>%
    mutate(
      method = method_name,
      # convert selected variable names into a logical indicator for all variables
      # TRUE = selected by method, FALSE = not selected -> used to compute TP,FP,etc
      selected = variable %in% selected_vars,
      beta_true = beta_true[variable]) %>%
    left_join(refit, by = "variable") %>%
    ## this coverage calculation is saying, "“If we fail to select a true 
    #  variable, it counts as a failure. Per Carter Sevick
    mutate(
      estimate_final = ifelse(selected, estimate, 0),
      bias = estimate_final - beta_true,
      covered = case_when(
        selected ~ (ci_low <= beta_true & ci_high >= beta_true),
        !selected & beta_true == 0 ~ TRUE,
        !selected & beta_true != 0 ~ FALSE),
      
    ## coverage is conditional on selection for the below code
    # mutate(
    #   bias = ifelse(selected, estimate - beta_true, NA_real_),
    #   covered = ifelse(
    #     selected,
    #     # **REVIEW** - not sure if this captures coverage correctly...
    #     ci_low <= beta_true & ci_high >= beta_true, NA),
      
      significant = ifelse(selected, p_value < alpha_test, FALSE),
    # β=0 AND selected AND significant
      type1_error = beta_true == 0 & selected & significant,
    # β≠0 AND not selected OR not significant
      # missed detection can happen by not selecting or finding significance
      type2_error = beta_true != 0 & (!selected | !significant))
  
  ######--- Confusion matrix ---######
  
  # number truly associated variables (beta ≠ 0) that were correctly selected
  TP <- sum(variable_level$selected & variable_level$beta_true != 0)
  # number of null variables (beta = 0) that were incorrectly selected
  FP <- sum(variable_level$selected & variable_level$beta_true == 0)
  # nmumber of truly associated variables (beta ≠ 0) that were NOT selected
  FN <- sum(!variable_level$selected & variable_level$beta_true != 0)
  # of null variables (beta = 0) that were correctly NOT selected
  TN <- sum(!variable_level$selected & variable_level$beta_true == 0)
  
  scenario_level <- tibble(
    method = method_name,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    false_positive_rate = FP / (FP + TN),
    false_negative_rate = FN / (FN + TP),
    model_size = length(selected_vars))
  
  list(
    variable_level = variable_level,
    scenario_level = scenario_level)
  }

mcse_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
  }

mcse_prop <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  p_hat <- mean(x)
  sqrt(p_hat * (1 - p_hat) / length(x))
  }


## Run Single Replication =================================================

run_single_rep <- function(n, rho, sim_id) {
  # simulates one dataset using hdrm
  dat <- hdrm::gen_data(
    n = n,
    p = p,
    p1 = p_true, # 5 true predictors
    beta = beta_true,
    family = "gaussian",
    corr = "exchangeable",
    rho = rho)
  
  # extract X and y from dat
  y <- dat$y
  X <- dat$X
  
  # store simulated data into df (optional, for backward selection)
  df <- data.frame(y = y, as.data.frame(X))
  # includes all 20 predictors
  full_model <- lm(y ~ ., data = df)
 
  ## Apply variable selection methods ##
  
  # backward selection using p-value / F-test 
  fit_p <- backward_pvalue(full_model, alpha_remove = alpha_pvalue)
  selected_p <- get_selected_lm(fit_p)
  
  # backward selection using AIC
  fit_aic <- MASS::stepAIC(full_model, direction = "backward", trace = FALSE)
  # chooses the model that minimizes AIC 
  selected_aic <- get_selected_lm(fit_aic)
  
  # ackward selection using BIC
  fit_bic <- MASS::stepAIC(
    full_model,
    direction = "backward",
    k = log(n),
    trace = FALSE)
  # chooses the model that minimizes BIC 
  selected_bic <- get_selected_lm(fit_bic)
  
  ######---Penalized Regression---######
  X_mat <- as.matrix(X)
  
  # LASSO
  cv_lasso <- cv.glmnet(
    x = X_mat,
    y = y,
    alpha = 1,
    nfolds = cv_folds,
    standardize = TRUE)
  
  # Elastic Net
  cv_enet <- cv.glmnet(
    x = X_mat,
    y = y,
    alpha = alpha_enet,
    nfolds = cv_folds,
    standardize = TRUE)
  
  selected_lasso_min <- get_selected_glmnet(cv_lasso, "lambda.min")
  selected_lasso_1se <- get_selected_glmnet(cv_lasso, "lambda.1se")
  
  selected_enet_min <- get_selected_glmnet(cv_enet, "lambda.min")
  selected_enet_1se <- get_selected_glmnet(cv_enet, "lambda.1se")
  
  # store selected variables from each method
  methods <- list(
    backward_pvalue = selected_p,
    backward_AIC = selected_aic,
    backward_BIC = selected_bic,
    lasso_lambda_min = selected_lasso_min,
    lasso_lambda_1se = selected_lasso_1se,
    elastic_net_lambda_min = selected_enet_min,
    elastic_net_lambda_1se = selected_enet_1se)
  
  ######---Evaluation of each method---######
  evals <- list()
  
  for (method_name in names(methods)) {
    selected_vars <- methods[[method_name]]
    evals[[method_name]] <- evaluate_method(
      method_name = method_name,
      selected_vars = selected_vars,
      y = y,
      X = X,
      beta_true = beta_true)
    }
  
  variable_results <- bind_rows(
    lapply(evals, function(x) x$variable_level)) %>%
    mutate(n = n, rho = rho, sim_id = sim_id)
  
  scenario_results <- bind_rows(
    lapply(evals, function(x) x$scenario_level)) %>%
    mutate(n = n, rho = rho, sim_id = sim_id)
  
  list(
    variable_results = variable_results, 
    scenario_results = scenario_results)
  }
  

## Main simulation loop =================================================

all_variable_results <- list()
all_scenario_results <- list()

counter <- 1

for (n in n_values) {
  for (rho in rho_values) {
    for (sim_id in 1:n_sim) {
      
      if (sim_id %% 50 == 0) {
        message("Running n = ", n, ", rho = ", rho, ", sim = ", sim_id)
        }
      
      out <- run_single_rep(n = n, rho = rho, sim_id = sim_id)
      
      all_variable_results[[counter]] <- out$variable_results
      all_scenario_results[[counter]] <- out$scenario_results
      
      counter <- counter + 1
    }
  }
  }

variable_results <- bind_rows(all_variable_results)
scenario_results <- bind_rows(all_scenario_results)


## Main simulation loop w/ Parallelization ====================================

t0 <- Sys.time()

handlers(global = TRUE)
handlers("progress") # or "txtprogressbar"

sim_grid <- expand.grid(
  n = n_values,
  rho = rho_values,
  sim_id = 1:n_sim)

# 18 cores but R crashes with 15 or more
plan(multisession, workers = parallel::detectCores() - 6)

with_progress({
  # progress bar
  prog <- progressor(along = seq_len(nrow(sim_grid)))
  # for each row of sim_grid, run one simulation replication.
  all_results <- future_lapply(
    seq_len(nrow(sim_grid)), # nrow(sim_grid)
    function(i) {
      prog()
      run_single_rep(
        n = as.numeric(sim_grid$n[i]), #as.numeric(sim_grid$n[i]),
        rho = as.numeric(sim_grid$rho[i]), # as.numeric(sim_grid$rho[i]),
        sim_id = as.integer(sim_grid$sim_id[i])) #as.integer(sim_grid$sim_id[i]))
      },
    future.seed = 2522) # review this******
  })

# simple run
with_progress({
  # progress bar
  prog <- progressor(along = seq_len(nrow(sim_grid)))
  # for each row of sim_grid, run one simulation replication.
  all_results <- future_lapply(
    seq_len(100), # nrow(sim_grid)
    function(i) {
      prog()
      run_single_rep(
        n = as.numeric(100), 
        rho = as.numeric(0.75), 
        sim_id = as.integer(1)) 
    },
    future.seed = 2522) 
})

plan(sequential)

variable_results <- bind_rows(
  lapply(all_results, function(x) x$variable_results))

scenario_results <- bind_rows(
  lapply(all_results, function(x) x$scenario_results))

# save - rather large file so use seriallizatin
qsave(
  all_results, 
  file.path(
    "/Users/ryan_summers/Library/Mobile\ Documents/com~apple~CloudDocs",
    "GitHub/BIOS6624-Class/Project4/DataProcessed/all_results.qs"), 
  preset = "fast", 
  nthreads = 8)

saveRDS(variable_results, "variable_results.rds")
saveRDS(scenario_results, "scenario_results.rds")

t1 <- Sys.time()

message("Runtime: ", 
        round(as.numeric(difftime(t1, t0, units = "mins")), 2), 
        " minutes")

# 17min w/12 scores on M5 max

## Summaries ====================================

# Method-level summary: selection accuracy and model size
method_summary <- scenario_results %>%
  group_by(n, rho, method) %>%
  summarise(
    # mean number of TPs across simulations
    mean_TP = mean(TP, na.rm = TRUE),
    # Monte Carlo SE (MCSE) of mean TP - uses variability across sim runs
    mcse_TP = mcse_mean(TP),
    # mean number of FPs
    mean_FP = mean(FP, na.rm = TRUE),
    # MCSE of mean FP
    mcse_FP = mcse_mean(FP),
    # avg TPR (sensitivity) across ims - prop of true vars correctly selected
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    # MCSE of sens
    mcse_sensitivity = mcse_mean(sensitivity),
    #avg specificity (TNR) across sims - prop of null vars correctly not selected 
    mean_specificity = mean(specificity, na.rm = TRUE),
    # MCSE of spec
    mcse_specificity = mcse_mean(specificity),
    # avg FPR across sims - prop of null vars incorrectly selected
    mean_false_positive_rate = mean(false_positive_rate, na.rm = TRUE),
    mcse_false_positive_rate = mcse_mean(false_positive_rate),
    # avg number of variables selected
    mean_model_size = mean(model_size, na.rm = TRUE),
    # MCSE of size
    mcse_model_size = mcse_mean(model_size),
    .groups = "drop")

# Variable-level summary: selection probability, bias, coverage, type I/II
variable_summary <- variable_results %>%
  group_by(n, rho, method, variable, beta_true) %>%
  summarise(
    # prob that a given variable is selected by the method across sim runs
    selection_probability = mean(selected),
    mcse_selection_probability = mcse_prop(selected),
    
    mean_bias = mean(bias, na.rm = TRUE),
    mcse_bias = mcse_mean(bias),
    
    mse = mean(bias^2, na.rm = TRUE),
    mcse_mse = mcse_mean(bias^2),
    
    coverage = mean(covered, na.rm = TRUE),
    mcse_coverage = mcse_prop(covered),
    # Type I error: only defined for the 15 null variables
    # Pr(reject H0 | beta = 0)
    type1_error_rate = ifelse(
      unique(beta_true) == 0,
      mean(type1_error), # median??
      NA_real_),
    mcse_type1_error = ifelse(
      unique(beta_true) == 0,
      mcse_prop(type1_error),
      NA_real_),
    
    # Type II error: only defined for the 5 true variables
    # Pr(fail to reject H0 | beta != 0)
    type2_error_rate = ifelse(
      unique(beta_true) != 0,
      mean(type2_error),
      NA_real_),
    mcse_type2_error = ifelse(
      unique(beta_true) != 0,
      mcse_prop(type2_error),
      NA_real_),
    .groups = "drop")

## Save to csv ====================================

# write.csv(method_summary, "method_summary.csv", row.names = F)
# write.csv(variable_summary, "variable_summary.csv", row.names = F)
# write.csv(scenario_results, "scenario_level_raw_results.csv", row.names = F)
# write.csv(variable_results, "variable_level_raw_results.csv", row.names = F)


## Sample Plots ====================================

# Metrics vs rho
ggplot(method_summary, aes(x = rho, y = mean_sensitivity, color = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ n) +
  scale_color_discrete(
    name = NULL,
    labels = c(
      backward_AIC = "AIC",
      backward_BIC = "BIC",
      backward_pvalue = "p-value",
      lasso_lambda_min = expression("LASSO (" * lambda[min] * ")"),
      lasso_lambda_1se = expression("LASSO (" * lambda[1*se] * ")"),
      elastic_net_lambda_min = expression("ElNet (" * lambda[min] * ")"),
      elastic_net_lambda_1se = expression("ElNet (" * lambda[1*se] * ")"))) +
  labs(y = "Sensitivity (TPR)", x = "Correlation (rho)") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(0.8, "lines"))

ggplot(variable_summary %>% filter(beta_true != 0),
       aes(x = method, y = mean_bias, color = method)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  facet_grid(n ~ rho) +
  scale_x_discrete(labels = c(
    backward_AIC = "AIC",
    backward_BIC = "BIC",
    backward_pvalue = "p-value",
    lasso_lambda_min = expression("LASSO (" * lambda[min] * ")"),
    lasso_lambda_1se = expression("LASSO (" * lambda[1*se] * ")"),
    elastic_net_lambda_min = expression("ElNet (" * lambda[min] * ")"),
    elastic_net_lambda_1se = expression("ElNet (" * lambda[1*se] * ")"))) +
  labs(
    title = "Bias of Est Coefficients (True Variables Only)",
    y = "Mean Bias",
    x = "Method") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))

# mean bias under Null and True vars
variable_results %>%
  mutate(var_type = ifelse(beta_true == 0, "Null", "True")) %>%
  group_by(var_type) %>%
  summarise(mean_bias = mean(bias, na.rm = TRUE) , .groups = "drop")

# Null bias relationship (in aggregate to just get a feel)
variable_results %>%
  filter(beta_true == 0) %>%
  summarise(
    total_null = n(),
    zero_bias_count = sum(bias == 0, na.rm = TRUE),
    prop_zero_bias = mean(bias == 0, na.rm = TRUE)) # proportion selected too


# bias across simulations by method
ggplot(variable_results %>% filter(beta_true != 0),
       aes(x = method, y = bias, fill = method)) +
  geom_boxplot() +
  facet_grid(n ~ rho) +
  scale_x_discrete(labels = c(
    backward_AIC = "AIC",
    backward_BIC = "BIC",
    backward_pvalue = "p-value",
    lasso_lambda_min = expression("LASSO (" * lambda[min] * ")"),
    lasso_lambda_1se = expression("LASSO (" * lambda[1*se] * ")"),
    elastic_net_lambda_min = expression("ElNet (" * lambda[min] * ")"),
    elastic_net_lambda_1se = expression("ElNet (" * lambda[1*se] * ")"))) +
  labs(
    title = "Distribution of Bias Across Simulations",
    y = "Bias",
    x = "Method") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))

# Type I error plotting
type1_summary <- variable_results %>%
  filter(beta_true == 0) %>%
  group_by(n, rho, method) %>%
  summarise(
    type1_error_rate = mean(type1_error, na.rm = TRUE), 
    mcse_type1_error = mcse_prop(type1_error),
    .groups = "drop")

ggplot(type1_summary,
       aes(x = method, y = type1_error_rate, fill = method)) +
  geom_col() +
  facet_grid(n ~ rho) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = c(
    backward_AIC = "AIC",
    backward_BIC = "BIC",
    backward_pvalue = "p-value",
    lasso_lambda_min = expression("LASSO (" * lambda[min] * ")"),
    lasso_lambda_1se = expression("LASSO (" * lambda[1*se] * ")"),
    elastic_net_lambda_min = expression("ElNet (" * lambda[min] * ")"),
    elastic_net_lambda_1se = expression("ElNet (" * lambda[1*se] * ")"))) +
  theme_bw() +
  labs(
    x = "Method",
    y = "Type I Error Rate",
    title = "Type I Error Rate Among Null Variables") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

# Type 1 by variable and method
ggplot(variable_summary %>% filter(beta_true == 0),
       aes(x = variable, y = type1_error_rate, color = method)) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(n ~ rho) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_color_manual(
    name = NULL,
    values = c(
      backward_AIC = "#E69F00",
      backward_BIC = "#56B4E9",
      backward_pvalue = "#009E73",
      elastic_net_lambda_min = "#0072B2",
      elastic_net_lambda_1se = "#CC79A7",
      lasso_lambda_min = "#D55E00",
      lasso_lambda_1se = "#999999"),
    labels = c(
      backward_AIC = "AIC",
      backward_BIC = "BIC",
      backward_pvalue = "p-value",
      lasso_lambda_min = expression("LASSO (" * lambda[min] * ")"),
      lasso_lambda_1se = expression("LASSO (" * lambda[1*se] * ")"),
      elastic_net_lambda_min = expression("ElNet (" * lambda[min] * ")"),
      elastic_net_lambda_1se = expression("ElNet (" * lambda[1*se] * ")"))) +
  labs(
    y = "Type I Error Rate",
    x = "Variable") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(variable_summary %>% filter(beta_true == 0),
       aes(x = method, y = type1_error_rate, fill = method)) +
  geom_boxplot() +
  facet_grid(n ~ rho) +
  scale_x_discrete(labels = c(
    backward_AIC = "AIC",
    backward_BIC = "BIC",
    backward_pvalue = "p-value",
    lasso_lambda_min = expression("LASSO (" * lambda[min] * ")"),
    lasso_lambda_1se = expression("LASSO (" * lambda[1*se] * ")"),
    elastic_net_lambda_min = expression("ElNet (" * lambda[min] * ")"),
    elastic_net_lambda_1se = expression("ElNet (" * lambda[1*se] * ")"))) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


ggplot(variable_summary %>% filter(beta_true == 0),
       aes(x = variable, y = type1_error_rate)) +
  geom_col() +
  facet_grid(method ~ n + rho) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Sample Tables ====================================

error_table <- variable_summary %>%
  group_by(n, rho, method) %>%
  summarise(
    `Type I Error (Null)` = mean(type1_error_rate[beta_true == 0], 
                                 na.rm = TRUE),
    `Type II Error (True)` = mean(type2_error_rate[beta_true != 0], 
                                  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    Method = recode(
      method,
      backward_pvalue = "Backward (p-value)",
      backward_AIC = "AIC",
      backward_BIC = "BIC",
      lasso_lambda_min = expression("LASSO (" * lambda[min] * ")"),
      lasso_lambda_1se = expression("LASSO (" * lambda[1*se] * ")"),
      elastic_net_lambda_min = expression("ElNet (" * lambda[1*se] * ")"),
      elastic_net_lambda_1se = expression("ElNet (" * lambda[1*se] * ")")
      )) %>%
  select(n, rho, Method, `Type I Error (Null)`, `Type II Error (True)`) %>% 
  gt() %>%
  tab_header(title = "Table: Type I & II Error Rates") %>%
  fmt_number(
    columns = c(`Type I Error (Null)`, `Type II Error (True)`),
    decimals = 3)







error_table_wide <- variable_summary %>%
  group_by(n, rho, method) %>%
  summarise(
    type1 = mean(type1_error_rate[beta_true == 0], na.rm = TRUE),
    type2 = mean(type2_error_rate[beta_true != 0], na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    Method = recode(
      method,
      backward_pvalue = "Backward (p-value)",
      backward_AIC = "AIC",
      backward_BIC = "BIC",
      lasso_lambda_min = "LASSO (λ<sub>min</sub>)", #(λₘᵢₙ) ; (λ₁ₛₑ)"
      lasso_lambda_1se = "LASSO (λ<sub>1se</sub>)",
      elastic_net_lambda_min = "ElNet (λ<sub>min</sub>))",
      elastic_net_lambda_1se = "ElNet (λ<sub>1se</sub>)"
    )) %>% 
  select(n, rho, Method, type1, type2) %>%
  pivot_wider(
    names_from = n,
    values_from = c(type1, type2),
    names_glue = "{.value}_{n}") %>%
  arrange(rho, Method)

error_table_wide %>%
  gt(rowname_col = "Method", groupname_col = "rho") %>%
  fmt_markdown(columns = Method) %>% 
  tab_header(title = "Type I & II Error Rates by Sample Size") %>%
  tab_spanner(
    label = "n = 250",
    columns = c(type1_250, type2_250)) %>%
  tab_spanner(
    label = "n = 500",
    columns = c(type1_500, type2_500)) %>%
  cols_label(
    type1_250 = "Type I",
    type2_250 = "Type II",
    type1_500 = "Type I",
    type2_500 = "Type II") %>%
  fmt_number(
    columns = starts_with("type"),
    decimals = 3)









