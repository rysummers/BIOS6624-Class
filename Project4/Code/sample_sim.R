
# Project 4 Simulation Skeleton ============================================
# Variable Selection Algorithms in Linear Regression

# Install if needed:
# install.packages(c("glmnet", "MASS", "dplyr", "purrr", "tibble"))
# remotes::install_github("pbreheny/hdrm")

library(hdrm)
library(glmnet)
library(MASS)
library(dplyr)
library(purrr)
library(tibble)
library(future)
library(furrr)
library(progress)


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

alpha_pvalue <- 0.15
alpha_enet <- 0.5
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
      selected = variable %in% selected_vars,
      beta_true = beta_true[variable]) %>%
    left_join(refit, by = "variable") %>%
    mutate(
      bias = ifelse(selected, estimate - beta_true, NA_real_),
      covered = ifelse(
        selected,
        ci_low <= beta_true & ci_high >= beta_true,
        NA),
      
      significant = ifelse(selected, p_value < alpha_test, FALSE),
      type1_error = beta_true == 0 & selected & significant,
      type2_error = beta_true != 0 & (!selected | !significant))
  
  TP <- sum(variable_level$selected & variable_level$beta_true != 0)
  FP <- sum(variable_level$selected & variable_level$beta_true == 0)
  FN <- sum(!variable_level$selected & variable_level$beta_true != 0)
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


## One simulation replicates =================================================

run_one_sim <- function(n, rho, sim_id) {
  # simulates one dataset using hdrm
  dat <- hdrm::gen_data(
    n = n,
    p = p,
    p1 = p_true, # 5 true predictors
    beta = beta_true,
    family = "gaussian",
    corr = "exchangeable",
    rho = rho)
  # store simulated data into df
  df <- data.frame(y = dat$y, as.data.frame(dat$X))
  # includes all 20 predictors
  full_model <- lm(y ~ ., data = df)
  
  # apply backward selection using p-value / F-test 
  fit_p <- backward_pvalue(full_model, alpha_remove = alpha_pvalue)
  selected_p <- get_selected_lm(fit_p)
  
  # apply backward selection using AIC
  fit_aic <- MASS::stepAIC(full_model, direction = "backward", trace = FALSE)
  selected_aic <- get_selected_lm(fit_aic)
  
  # apply backward selection using BIC
  fit_bic <- MASS::stepAIC(
    full_model,
    direction = "backward",
    k = log(n),
    trace = FALSE)
  # chooses the model that minimizes BIC 
  selected_bic <- get_selected_lm(fit_bic)
  
  ## apply Lasso and elastic net
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
  
  methods <- list(
    backward_pvalue = selected_p,
    backward_AIC = selected_aic,
    backward_BIC = selected_bic,
    lasso_lambda_min = selected_lasso_min,
    lasso_lambda_1se = selected_lasso_1se,
    elastic_net_lambda_min = selected_enet_min,
    elastic_net_lambda_1se = selected_enet_1se)
  
  evals <- imap(methods, ~ evaluate_method(.y, .x, y, X, beta_true))
  
  variable_results <- map_dfr(evals, "variable_level") %>%
    mutate(n = n, rho = rho, sim_id = sim_id)
  
  scenario_results <- map_dfr(evals, "scenario_level") %>%
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
      
      out <- run_one_sim(n = n, rho = rho, sim_id = sim_id)
      
      all_variable_results[[counter]] <- out$variable_results
      all_scenario_results[[counter]] <- out$scenario_results
      
      counter <- counter + 1
    }
  }
  }

variable_results <- bind_rows(all_variable_results)
scenario_results <- bind_rows(all_scenario_results)




## Main simulation loop w/ Parallelization ====================================

handlers(global = TRUE)
handlers("progress")   # or "txtprogressbar"

sim_grid <- expand.grid(
  n = n_values,
  rho = rho_values,
  sim_id = 1:n_sim)

plan(multisession, workers = parallel::detectCores() - 3)

with_progress({
  prog <- progressor(along = 1:nrow(sim_grid))
  all_results <- future_pmap(
    sim_grid,
    function(n, rho, sim_id) {
      prog()
      run_one_sim(
        n = as.numeric(n),
        rho = as.numeric(rho),
        sim_id = as.integer(sim_id))
      },
    .options = furrr_options(seed = TRUE))
  })

plan(sequential)

variable_results <- bind_rows(lapply(all_results, `[[`, "variable_results"))
scenario_results <- bind_rows(lapply(all_results, `[[`, "scenario_results"))

## Summaries ====================================

# Method-level summary: selection accuracy and model size
method_summary <- scenario_results %>%
  group_by(n, rho, method) %>%
  summarise(
    mean_TP = mean(TP),
    mcse_TP = mcse_mean(TP),
    
    mean_FP = mean(FP),
    mcse_FP = mcse_mean(FP),
    
    sensitivity = mean(sensitivity),
    mcse_sensitivity = mcse_mean(sensitivity),
    
    specificity = mean(specificity),
    mcse_specificity = mcse_mean(specificity),
    
    false_positive_rate = mean(false_positive_rate),
    mcse_false_positive_rate = mcse_mean(false_positive_rate),
    
    mean_model_size = mean(model_size),
    mcse_model_size = mcse_mean(model_size),
    
    .groups = "drop")

# Variable-level summary: selection probability, bias, coverage, type I/II
variable_summary <- variable_results %>%
  group_by(n, rho, method, variable, beta_true) %>%
  summarise(
    selection_probability = mean(selected),
    mcse_selection_probability = mcse_prop(selected),
    
    mean_bias_among_selected = mean(bias, na.rm = TRUE),
    mcse_bias_among_selected = mcse_mean(bias),
    
    mse_among_selected = mean(bias^2, na.rm = TRUE),
    mcse_mse_among_selected = mcse_mean(bias^2),
    
    coverage_among_selected = mean(covered, na.rm = TRUE),
    mcse_coverage_among_selected = mcse_prop(covered),
    
    type1_error_rate = ifelse(
      unique(beta_true) == 0,
      mean(type1_error),
      NA_real_),
    
    mcse_type1_error = ifelse(
      unique(beta_true) == 0,
      mcse_prop(type1_error),
      NA_real_),
    
    type2_error_rate = ifelse(
      unique(beta_true) != 0,
      mean(type2_error),
      NA_real_),
    
    mcse_type2_error = ifelse(
      unique(beta_true) != 0,
      mcse_prop(type2_error),
      NA_real_), .groups = "drop")

## Save results ====================================

write.csv(method_summary, "method_summary.csv", row.names = FALSE)
write.csv(variable_summary, "variable_summary.csv", row.names = FALSE)
write.csv(scenario_results, "scenario_level_raw_results.csv", row.names = FALSE)
write.csv(variable_results, "variable_level_raw_results.csv", row.names = FALSE)

