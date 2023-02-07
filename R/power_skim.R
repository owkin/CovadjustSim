#' Skim through the power range of a simulation
#'
#' Loosely compute the statistical power of the treatment effect at a series of
#' `n_samples` points, until the upper limit of power is reached. At each step,
#' the increment of `n_samples` is estimated adaptively.
#' @param x An object of a subclass of [simu_ph_base]
#' @param f A formula specifying the covariates
#' @param n_reps Number of repetitions of the trial simulation at each
#'   `n_samples` to estimate the statistical power of the treatment effect.
#'   Intended to be small to keep the speed of computation.
#' @param n_points Minimal number of `n_samples` points to be examined
#' @param n_init starting point of `n_samples`
#' @param power_stop Upper limit of power to stop the skim
#' @param ... further arguments passed to other methods.
#' @return A data frame of `n_samples`, `n_events` and `power`
#' @seealso [interpolate_power()]
#' @export
power_skim <- function(x, ...) {
  UseMethod("power_skim")
}

#' @rdname power_skim
#' @export
power_skim.simu_ph_base <- function(x, f, n_reps = 300, n_points = 7,
                                    n_init = 20, power_stop = 0.9, ...) {
  step <- 20  # a manually set conservative initial value
  step_max <- 200
  n_samples <- n_init - step
  results <- data.frame()
  repeat {
    n_samples <- n_samples + min(step_max, step)
    results <- results %>% bind_rows(x %>% compute_power(
      n_reps, n_samples, list(f), get_cindex = FALSE, get_lambda = FALSE
    ))
    k <- nrow(results)
    power_k <- results$power[k]
    # Stop condition
    if ((k >= n_points && power_k >= power_stop) || power_k > 0.95) {
      break
    }
    if (k == 1 || power_k >= power_stop) {
      next
    }
    # Update step
    ratio <- (power_k - results$power[1]) /
        (results$n_samples[k] - results$n_samples[1])
    if (ratio > 0) {
      n_remained <- max(1, n_points - k)
      step <- ceiling((power_stop - power_k) / (ratio * n_remained))
    }
  }
  return(results)
}

#' Find values of `n_samples` for given values of statistical power
#'
#' Often combined with the function [power_skim()] to get a quick estimation of
#' the required sample size and number of events for a treatment effect
#' estimation to reach a predefined power.
#'
#' @param df_power A data frame of `power`, `n_samples` and `n_events`
#' @param power_targets Values of `power` for which corresponding `n_samples`
#'   and `n_events` will interpolated
#' @seealso [power_skim()]
#' @export
interpolate_power <- function(df_power, power_targets) {
  model_sample <- lm(power ~ poly(n_samples, 2), data = df_power)
  model_event <- lm(power ~ poly(n_events, 2), data = df_power)
  df_sample <- data.frame(n_samples = 1:max(df_power$n_samples)) %>%
    mutate(power = predict(model_sample, .))
  df_event <- data.frame(n_events = 1:max(df_power$n_events)) %>%
    mutate(power = predict(model_event, .))
  indices_sample <- sapply(power_targets, function(power) {
    return(which.min(abs(df_sample$power - power)))
  })
  indices_event <- sapply(power_targets, function(power) {
    return(which.min(abs(df_event$power - power)))
  })
  return(data.frame(
    power = power_targets,
    n_samples = df_sample$n_samples[indices_sample],
    n_events = df_event$n_events[indices_event]
  ))
}
