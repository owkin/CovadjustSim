#' Estimate the reduction in sample size required to reach 80% statistical power
#' when adjusting for a covariate in a clinical trial.
#'
#' For each formula, a list of `n_samples` will be generated. For each
#' `n_samples`, statistical power of the treatment effect is estimated. Then for
#' each formula, the `power ~ n_samples` curve will be interpolated to estimate
#' the sample size required to achieve 80% power.
#'
#' @param simu An object of class [simu_ph_weibull]
#' @param formula1 The baseline formula to be analysed in `simu`
#' @param formula2 The formula to be compared with `formula1` to estimate the
#'   reduction in sample size to reach 80% power
#' @param begin1,end1,begin2,end2 Values of `n_samples` specifying the range of
#'   the list of `n_samples` for each formula, often fixed by [power_skim()].
#' @param n_reps Number of repeated experiments to estimate the statistical
#'   power of the treatment effect.
#' @param n_points Minimum number of sample size points used to interpolate the
#'   sample size that gives statistical power of 80%.
#' @param seed A random seed for reproducibility.
#' @param n_cores Number of cores to be used in parallel computing.
#' @return A data frame including the average cumulative incidence `lambda`, the
#'   average concordance `cindex` for each formula, the number of samples
#'   corresponding to 80% statistical power for each formula, and the estimated
#'   sample size reduction.
#' @examples
#' simu <- simu_ph_weibull(
#'   shape = 0.5, hr_trt = 0.4, dropout = 0.01, t_cutoff = 5, beta = 0.58,
#'   kappa = -1.53
#' )
#' formula1 <- survival::Surv(time, event) ~ 1
#' formula2 <- survival::Surv(time, event) ~ X
#'
#' estimate_sample_size_reduction(simu,
#'   formula1, 80, 200,
#'   formula2, 70, 170,
#'   n_reps = 30, n_points = 5, seed = 1
#' )
#' @seealso [power_skim()]
#' @export
estimate_sample_size_reduction <- function(simu, formula1, begin1, end1,
                                           formula2, begin2, end2,
                                           n_reps = 1e+3, n_points = 6,
                                           seed = NULL, n_cores = 1) {
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
  df_power1 <- purrr::map_dfr(ceiling(seq(begin1, end1, length.out = n_points)),
    function(n_samples) {
      return(simu %>% compute_power(n_reps, n_samples, list(formula1)))
    }
  )
  n80_1 <- interpolate_power(df_power1, 0.8)$n_samples
  df_power2 <- purrr::map_dfr(ceiling(seq(begin2, end2, length.out = n_points)),
    function(n_samples) {
      return(simu %>% compute_power(n_reps, n_samples, list(formula2)))
    }
  )
  n80_2 <- interpolate_power(df_power2, 0.8)$n_samples

  results <- data.frame(
    lambda = mean(c(df_power1$lambda, df_power2$lambda), na.rm = TRUE),
    cindex1 = mean(df_power1$cindex, na.rm = TRUE),
    cindex2 = mean(df_power2$cindex, na.rm = TRUE),
    n80_1 = n80_1,
    n80_2 = n80_2,
    r2_obs_80 = 1 - n80_2 / n80_1
  )
  return(results)
}
