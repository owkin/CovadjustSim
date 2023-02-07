#' Generate a dataset simulating a clinical trial
#'
#' No default method. Each class (subclass of [simu_ph_base]) should implement
#' this method.
#' @param x An object of a subclass of [simu_ph_base].
#' @param n_samples Number of samples in the returned data frame.
#' @param ... further arguments passed to other methods.
#' @return A data frame of simulated clinical trial with at least three columns:
#'   time: numeric, the survival time.
#'   event: integer 1 or 0, indicating either an event (e.g. death) or right
#'     censoring.
#'   treatment: integer 1 or 0, indicating whether the patient is in the
#'     treatment arm or in the control arm.
#'   Additionally, other columns can exist which are considered as potential
#'   covariates in the study.
#' @export
generate_trial_data <- function(x, ...) {
  UseMethod("generate_trial_data")
}

#' @rdname generate_trial_data
#' @export
generate_trial_data.simu_ph_weibull <- function(x, n_samples, ...) {
  covariate <- rnorm(ceiling(n_samples / x$inclusion))
  threshold <- quantile(covariate, x$inclusion)
  covariate <- covariate[covariate <= threshold] %>% head(n = n_samples)

  treatment <- rep(sample(0:1), length.out = n_samples)
  cox_linear <- x$beta * covariate + x$kappa + log(x$hr_trt) * treatment
  scale <- exp(-cox_linear)^(1 / x$shape)

  data <- generate_survival_data(
    n_samples, x$shape, scale, x$dropout, x$t_cutoff
  ) %>%
    mutate(
      X = covariate,
      treatment = treatment
    )
  return(data)
}

generate_survival_data <- function(n_samples, shape, scale, dropout, t_cutoff) {
  time_censored <- NULL
  if (dropout > 0) {
    time_censored <- pmin(t_cutoff, rexp(n_samples, dropout))
  } else {
    time_censored <- rep(t_cutoff, n_samples)
  }
  time_event <- rweibull(n_samples, shape, scale)

  return(data.frame(
    time = pmin(time_event, time_censored),
    event = as.integer(time_event <= time_censored)
  ))
}
