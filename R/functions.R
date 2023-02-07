#' Fix `kappa` and `beta` parameters for multiple experiment settings
#'
#' The relations between (`cindex`, `lambda`) and (`kappa`, `beta`) are
#' determined given a tuple (`shape`, `dropout`, `t_cutoff`).
#' For each of these tuples, the corresponding (`kappa`, `beta`) pairs for
#' multiple (`cindex`, `lambda`) pairs can be interpolated with the same model.
#' First, a fine grid of `beta` and `kappa` (ranges depend on `shape` and the
#' covariate, and are fixed manually) is generated, the values of each (cindex,
#' lambda) pair is estimated for each point on the grid. Then the value of
#' `beta` corresponding to the desired `cindex` is interpolated (`cindex` does
#' not depend on `kappa`) with a polynomial. Then given the values of `beta` and
#' `lambda`, the value of `kappa` is interpolated with [interp::interp()].
#'
#' @param param_grid A data frame. Must contains the following columns:
#'   `shape`, `dropout`, `t_cutoff`, `cindex`, `lambda`.
#' @param ... Additional arguments passed on to `get_mapping_grid`.
#' @param seed A random seed for reproducibility.
#' @param n_cores Number of cores to be used in parallel computing.
#' @return A data frame of `param_grid` with two additional columns: `beta` and
#'   `kappa`.
#' @examples
#' param_grid <- tidyr::expand_grid(
#'   shape = 0.5,
#'   dropout = 0.01,
#'   t_cutoff = 5,
#'   n_reps = 10000,
#'   hr_trt = 0.7,
#'   cindex = c(0.75, 0.85),
#'   lambda = seq(0.1, 0.6, 0.1)
#' )
#' beta_list <- c(0.001, seq(0.1, 3, length.out = 9))
#' kappa_list <- seq(-8 + 0.4, 4, length.out = 7)
#' fix_dist_params(
#'   param_grid, beta_list = beta_list, kappa_list = kappa_list,
#'   n_samples = 500, seed = 42
#' )
#' @export
fix_dist_params <- function(param_grid, ..., seed = NULL, n_cores = 1) {
  params <- list(...)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (n_cores > 1) {
    strategy <- multisession
    if (parallelly::supportsMulticore()) {
      strategy <- multicore
    }
    withr::local_(plan)(strategy, workers = n_cores)
  }
  res <- param_grid %>%
    group_by(shape, dropout, t_cutoff) %>%
    group_modify(~ {
      # .y: group keys, .x: grouped data frame
      mapping <- do.call(get_mapping_grid, c(.y, params))
      df_mapped <- purrr::map2_dfr(
        .x$cindex, .x$lambda, interpolate_params, mapping
      )
      return(.x %>% bind_cols(df_mapped))
    }) %>%
    ungroup()
  return(res)
}

#' Find by interpolation `beta` and `kappa` that generate dataset with desired
#' `cindex` and `lambda`
#'
#' `cindex` is assumed to be dependent mainly on `beta` and slightly on `kappa`.
#' Interpolation is done based on the dataframe `mapping`, and by iteration.
#' `mapping` has four columns: `beta`, `kappa`, `cindex`, `lambda`.
#' In each iteration, `beta_pred` is first interpolated by a `beta ~ cindex`
#' curve at `cindex_target`, `kappa_pred` is then interpolated by the
#' `kappa ~ beta, lambda` plane at (`beta_pred`, `lambda_target').
#' @noRd
interpolate_params <- function(cindex_target, lambda_target, mapping) {
  mapping <- as.data.frame(mapping) %>% drop_na()

  # Initialize `beta ~ cindex` curve by averaging over `kappa`
  df_beta <- mapping %>%
    group_by(beta) %>%
    summarise(cindex = mean(cindex))

  # Number of iterations fixed manually
  for (i in 1:3) {
    beta_pred <- lm(beta ~ poly(cindex, degree = 5), data = df_beta) %>%
      predict(data.frame(cindex = cindex_target))
    kappa_pred <- interp::interp(
      mapping$beta, mapping$lambda, mapping$kappa, beta_pred, lambda_target,
      duplicate = "user", dupfun = function(x) x[[which.min(abs(x))]],
      output = "points"
    )$z
    # Update `beta ~ cindex` curve by slicing at `kappa` closest to `kappa_pred`
    df_beta <- mapping %>%
      filter(abs(kappa - kappa_pred) == min(abs(kappa - kappa_pred)))
  }

  if (is.null(beta_pred) || is.null(kappa_pred)) {
    warning(paste(
      "Search of `kappa` and `beta` failed at cindex", round(cindex_target, 3),
      "and lambda", round(lambda_target, 3), "failed."
    ))
  }

  return(data.frame(beta = beta_pred, kappa = kappa_pred))
}

get_mapping_grid <- function(shape, dropout, t_cutoff,
                             beta_list = c(0.001, seq(0.1, 3, 0.1)),
                             kappa_list = seq(-8 + 0.4, 4, 0.4),
                             n_samples = 5000, seed = NULL, n_cores = 1) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (n_cores > 1) {
    strategy <- multisession
    if (parallelly::supportsMulticore()) {
      strategy <- multicore
    }
    withr::local_(plan)(strategy, workers = n_cores)
  }

  param_grid <- tidyr::expand_grid(beta = beta_list, kappa = kappa_list)
  res <- future_pmap_dfr(
    param_grid, get_param_mapping, shape, dropout, t_cutoff, n_samples,
    .options = furrr_options(seed = TRUE)
  )
  # To avoid duplicated points during [interp::interp()] due to precision issue,
  # where two points are considered distinct in R, but duplicated in C++.
  res <- res %>% mutate(across(everything(), ~ round(.x, 6)))
  return(res)
}

get_param_mapping <- function(beta, kappa, shape, dropout, t_cutoff,
                              n_samples) {
  res <- data.frame(kappa = kappa, beta = beta, lambda = NA, cindex = NA)

  # Set hr_trt = 1 to return only the control arm
  simu <- simu_ph_weibull(shape, hr_trt = 1, dropout, t_cutoff, beta, kappa)
  data <- simu %>% generate_trial_data(n_samples = n_samples)

  if (sum(data$event) == 0) {
    return(res)
  }

  km <- survfit(Surv(time, event) ~ 1, data = data)
  res$lambda <- 1 - tail(km$surv, 1)
  res$cindex <- coxph(Surv(time, event) ~ X, data = data)$concordance[[6]]

  return(res)
}
