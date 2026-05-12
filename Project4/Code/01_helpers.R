# Helper functions for Project 4 simulation

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
        data = model.frame(current_model)
      )
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
    return(tibble::tibble(
      variable = character(),
      estimate = numeric(),
      p_value = numeric(),
      ci_low = numeric(),
      ci_high = numeric()
    ))
  }
  
  form <- as.formula(paste("y ~", paste(selected_vars, collapse = " + ")))
  fit <- lm(form, data = df)
  
  coefs <- summary(fit)$coefficients
  cis <- confint(fit)
  
  tibble::tibble(
    variable = rownames(coefs)[-1],
    estimate = coefs[-1, "Estimate"],
    p_value = coefs[-1, "Pr(>|t|)"],
    ci_low = cis[-1, 1],
    ci_high = cis[-1, 2])
  }

evaluate_method <- function(method_name, selected_vars, y, X, beta_true, alpha_test = 0.05) {
  all_vars <- names(beta_true)
  refit <- post_selection_refit(y, X, selected_vars)
  
  variable_level <- tibble::tibble(variable = all_vars) %>%
    dplyr::mutate(
      method = method_name,
      selected = variable %in% selected_vars,
      beta_true = beta_true[variable]) %>%
    dplyr::left_join(refit, by = "variable") %>%
    dplyr::mutate(
      estimate_final = ifelse(selected, estimate, 0),
      bias = estimate_final - beta_true,
      covered = dplyr::case_when(
        selected ~ (ci_low <= beta_true & ci_high >= beta_true),
        !selected & beta_true == 0 ~ TRUE,
        !selected & beta_true != 0 ~ FALSE),
      significant = ifelse(selected, p_value < alpha_test, FALSE),
      type1_error = beta_true == 0 & selected & significant,
      type2_error = beta_true != 0 & (!selected | !significant))
  
  TP <- sum(variable_level$selected & variable_level$beta_true != 0)
  FP <- sum(variable_level$selected & variable_level$beta_true == 0)
  FN <- sum(!variable_level$selected & variable_level$beta_true != 0)
  TN <- sum(!variable_level$selected & variable_level$beta_true == 0)
  
  scenario_level <- tibble::tibble(
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
  
  list(variable_level = variable_level, scenario_level = scenario_level)
}

mcse_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

mcse_prop <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  p_hat <- mean(x)
  sqrt(p_hat * (1 - p_hat) / length(x))
}

mcse_bin <- function(p, nsim) {
  sqrt(p * (1 - p) / nsim)
}
