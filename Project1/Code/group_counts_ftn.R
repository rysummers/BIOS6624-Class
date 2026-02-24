
# function for calculating obs per grouping variable
count_grp_obs <- function(data, vars, include_na = TRUE, id = NULL) {
  stopifnot(is.data.frame(data))
  
  dat <- dplyr::as_tibble(data) %>%
    dplyr::select(dplyr::all_of(vars))
  
  dat_long <- dat %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "variable",
      values_to = "category") %>%
    dplyr::mutate(
      category = dplyr::if_else(is.na(category), 
                                "Missing", 
                                as.character(category)))
  
  if (!include_na) dat_long <- dplyr::filter(dat_long, category != "Missing")
  
  if (is.null(id)) {
    dat_long %>%
      dplyr::count(variable, category, name = "n") %>%
      dplyr::group_by(variable) %>%
      dplyr::mutate(pct = n / sum(n)) %>%
      dplyr::ungroup()
  } else {
    data %>%
      dplyr::select({{ id }}, dplyr::all_of(vars)) %>%
      tidyr::pivot_longer(-{{ id }}, 
                          names_to = "variable", 
                          values_to = "category") %>%
      dplyr::mutate(
        category = dplyr::if_else(is.na(category), 
                                  "Missing", 
                                  as.character(category))) %>%
      dplyr::distinct({{ id }}, variable, category) %>%
      dplyr::count(variable, category, name = "n_subjects") %>%
      dplyr::group_by(variable) %>%
      dplyr::mutate(pct_subjects = n_subjects / sum(n_subjects)) %>%
      dplyr::ungroup()
  }}