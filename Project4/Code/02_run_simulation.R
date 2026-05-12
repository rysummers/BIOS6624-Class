#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdrm)
  library(glmnet)
  library(MASS)
  library(dplyr)
  library(tibble)
  library(future)
  library(future.apply)
  library(progressr)
  library(here)
})

source(here("Code", "00_config.R"))
source(here("Code", "01_helpers.R"))

run_single_rep <- function(n, rho, sim_id) {
  dat <- hdrm::gen_data(
    n = n,
    p = p,
    p1 = p_true,
    beta = beta_true,
    family = "gaussian",
    corr = "exchangeable",
    rho = rho)
  
  y <- dat$y
  X <- dat$X
  df <- data.frame(y = y, as.data.frame(X))
  full_model <- lm(y ~ ., data = df)
  
  fit_p <- backward_pvalue(full_model, alpha_remove = alpha_pvalue)
  selected_p <- get_selected_lm(fit_p)
  
  fit_aic <- MASS::stepAIC(full_model, direction = "backward", trace = FALSE)
  selected_aic <- get_selected_lm(fit_aic)
  
  fit_bic <- MASS::stepAIC(full_model, direction = "backward", k = log(n), 
                           trace = FALSE)
  selected_bic <- get_selected_lm(fit_bic)
  
  X_mat <- as.matrix(X)
  
  cv_lasso <- cv.glmnet(
    x = X_mat,
    y = y,
    alpha = 1,
    nfolds = cv_folds,
    standardize = TRUE)
  
  cv_enet <- cv.glmnet(
    x = X_mat,
    y = y,
    alpha = alpha_enet,
    nfolds = cv_folds,
    standardize = TRUE)
  
  methods <- list(
    backward_pvalue = selected_p,
    backward_AIC = selected_aic,
    backward_BIC = selected_bic,
    lasso_lambda_min = get_selected_glmnet(cv_lasso, "lambda.min"),
    lasso_lambda_1se = get_selected_glmnet(cv_lasso, "lambda.1se"),
    elastic_net_lambda_min = get_selected_glmnet(cv_enet, "lambda.min"),
    elastic_net_lambda_1se = get_selected_glmnet(cv_enet, "lambda.1se"))
  
  evals <- vector("list", length(methods))
  names(evals) <- names(methods)
  
  for (method_name in names(methods)) {
    evals[[method_name]] <- evaluate_method(
      method_name = method_name,
      selected_vars = methods[[method_name]],
      y = y,
      X = X,
      beta_true = beta_true)
  }
  
  variable_results <- dplyr::bind_rows(lapply(evals, function(x) x$variable_level)) %>%
    dplyr::mutate(n = n, rho = rho, sim_id = sim_id)
  
  scenario_results <- dplyr::bind_rows(lapply(evals, function(x) x$scenario_level)) %>%
    dplyr::mutate(n = n, rho = rho, sim_id = sim_id)
  
  list(variable_results = variable_results, scenario_results = scenario_results)
}

message("Starting Project 4 simulation...")
message("n_sim = ", n_sim, "; alpha_pvalue = ", alpha_pvalue)
message("Worst-case MCSE for proportion metrics: ", round(mcse_bin(0.50, n_sim), 5))
message("MCSE near 95% coverage: ", round(mcse_bin(0.95, n_sim), 5))

t0 <- Sys.time()

sim_grid <- expand.grid(
  n = n_values,
  rho = rho_values,
  sim_id = seq_len(n_sim))

# for mac
#workers <- max(1, parallel::detectCores() - 3)
# for alpine
workers <- 18

message("Using workers: ", workers)
future::plan(future::multisession, workers = workers)
progressr::handlers(global = TRUE)
progressr::handlers("progress")

all_results <- progressr::with_progress({
  prog <- progressr::progressor(along = seq_len(nrow(sim_grid)))
  future.apply::future_lapply(
    seq_len(nrow(sim_grid)),
    function(i) {
      prog()
      run_single_rep(
        n = as.numeric(sim_grid$n[i]),
        rho = as.numeric(sim_grid$rho[i]),
        sim_id = as.integer(sim_grid$sim_id[i]))
    },
    future.seed = seed_value)
})

future::plan(future::sequential)

variable_results <- dplyr::bind_rows(
  lapply(all_results, function(x) x$variable_results))
scenario_results <- dplyr::bind_rows(
  lapply(all_results, function(x) x$scenario_results))

saveRDS(variable_results, here("DataProcessed","variable_results.rds"))
saveRDS(scenario_results, here("DataProcessed", "scenario_results.rds"))

# Optional: save full results only if needed; it is large.
# saveRDS(all_results, file.path(out_dir, "all_results.rds"))

t1 <- Sys.time()
message("Runtime: ", 
        round(as.numeric(difftime(t1, t0, units = "mins")), 2), 
        " minutes")
message("Saved results to: ", file.path(project_dir, "DataProcessed"))
