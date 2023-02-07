#' Generic S3 class of basic trial simulation based on proportional hazard model
#'
#' The class has no member and a method [compute_power()]. Derived classes
#' should implement the method [generate_trial_data()], and optionally the
#' method [compute_power()].
#'
#' @param class A string specifying the subclass
#' @return An object of class `simu_ph_base` or the derived subclass
#' @seealso [simu_ph_weibull]
#' @keywords internal
simu_ph_base <- function(class = character()) {
  simu <- list()
  class(simu) <- c(class, "simu_ph_base")
  return(simu)
}

#' S3 class of trial simulation based on Weibull and proportional hazard model
#'
#' Baseline survival times are modeled by the Weibull distribution.
#' Covariate is a variable `X` following the standard normal distribution.
#' Overall survival is modeled by the proportional hazard model.
#' Right censoring time is modeled by an exponential distribution with specified
#' dropout rate.
#' @param shape The shape parameter of the Weibull distribution.
#' @param hr_trt The treatment hazard ratio.
#' @param dropout The rate parameter in the exponential distribution to model
#'   the right censoring.
#' @param t_cutoff The cut-off time (in year) beyond which all patients will be
#'   right censored.
#' @param beta,kappa Parameters of the hazard term `kappa + beta * X` where
#'   `X` is the covariate. Their values are often fixed so as to reach a desired
#'   cumulative incidence `lambda` and a desired concordance `cindex' of the
#'   fitted proportional hazard model at `t_cutoff` in the control arm.
#' @param inclusion A numerical value between ]0, 1] to simulate the inclusion
#'   criterion. Specifically, it is a quantile threshold of the covariate `X`,
#'   only samples with a covariate value below the threshold are included in the
#'   trial.
#' @return An object of class `c(simu_ph_weibull, simu_ph_base)`
#' @seealso [simu_ph_base]
#' @export
simu_ph_weibull <- function(shape, hr_trt, dropout, t_cutoff, beta, kappa,
                       inclusion = 1) {
  object <- simu_ph_base(class = "simu_ph_weibull")
  object$shape <- shape
  object$hr_trt <- hr_trt
  object$dropout <- dropout
  object$t_cutoff <- t_cutoff
  object$kappa <- kappa
  object$beta <- beta
  object$inclusion <- inclusion

  return(object)
}
