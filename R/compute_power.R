#' Compute the statistical power of a treatment effect
#'
#' Compute the statistical power of a treatment effect from a simulated dataset
#' for a set of formulae, compute optionally also the concordance and cumulative
#' incidence.
#' @param x An object of a subclass of [simu_ph_base].
#' @param n_reps Number of repetitions of the trial simulation to estimate the
#'   statistical power of the treatment effect.
#' @param n_samples Number of samples in each trial simulation.
#' @param formulae A (named) list of formulae specifying the group of covariates
#'   used in the estimation of treatment effect. Each formula will be passed to
#'   [survival::coxph()], and should be of the form `Surv(time, event) ~
#'   covariate_1 + covariate_2 + ...`.
#' @param get_cindex Logical. If TRUE, estimate the concordance for each of the
#'   formula in `formulae` in the control arm (treatment == 0) of the trial.
#' @param get_lambda Logical. If TRUE, estimate the cumulative incidence of the
#'   event at the end of the trial in the control arm.
#' @param seed A random seed for reproducibility.
#' @param n_cores Number of cores to be used in parallel computing.
#' @param ... further arguments passed to other methods.
#' @return A data frame
#' @examples
#' simu <- simu_ph_weibull(
#'   shape = 0.5, hr_trt = 0.4, dropout = 0.1, t_cutoff = 5,
#'   beta = 1, kappa = -1, inclusion = 0.8
#' )
#' formulae <- list(
#'   base = survival::Surv(time, event) ~ 1,
#'   adj = survival::Surv(time, event) ~ X
#' )
#' res <- compute_power(simu, n_reps = 10, n_samples = 500, formulae = formulae)
#' @export
compute_power <- function(x, ...) {
  UseMethod("compute_power")
}

#' @rdname compute_power
#' @export
compute_power.simu_ph_base <- function(x, n_reps, n_samples, formulae,
                                    get_cindex = TRUE, get_lambda = TRUE,
                                    seed = NULL, n_cores = 1, ...) {
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

  boot <- future_map_dfr(1:n_reps,
    function(i) {
      return(x %>%
        generate_trial_data(n_samples) %>%
        single_trial(formulae, get_cindex, get_lambda)
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  return(boot %>%
    group_by(group) %>%
    summarise(
      n_samples  = n_samples,
      across(!"pvalue", ~ mean(.x, na.rm = TRUE)),
      across("pvalue", ~ mean(.x < 0.05, na.rm = TRUE))
    ) %>%
    rename("power" = "pvalue")
  )
}

single_trial <- function(data, formulae, get_cindex = TRUE, get_lambda = TRUE) {
  res <- data.frame(matrix(NA, nrow = length(formulae), ncol = 5))
  colnames(res) <- c("n_events", "cindex", "lambda", "pvalue", "group")

  if (!is.null(names(formulae))) {
    res["group"] <- factor(names(formulae))
  }
  n_events <- sum(data$event)
  res["n_events"] <- n_events
  if (n_events == 0) {
    return(res)
  }
  res["pvalue"] <- sapply(formulae, function(f) {
    tryCatch(
      {
        model <- coxph(update.formula(f, ~ treatment + .), data = data)
        # 5th column is `Pr(>|z|)`
        return(summary(model)$coefficients["treatment", 5])
      },
      error = function(e) {
        return(NA)
      }
    )
  })

  data_control <- NULL
  if (get_lambda || get_cindex) {
    data_control <- data %>% filter(treatment == 0)
  }

  if (get_lambda) {
    km <- survfit(Surv(time, event) ~ 1, data = data_control)
    res["lambda"] <- 1 - tail(km$surv, 1)
  }

  if (get_cindex) {
    res["cindex"] <- sapply(formulae, function(f) {
      tryCatch(
        {
          model <- coxph(f, data = data_control)
          return(model$concordance[["concordance"]])
        },
        error = function(e) {
          return(NA)
        }
      )
    })
  }

  return(res)
}
