#' Compute various measures of R-squared of the proportional hazard model for a
#' given dataset and group of covariates
#'
#' @param data A data frame of survival data
#' @param f A formula describing the group of covaraites
#' @param cs_only Logical. If TRUE, only the Cox-Snell R-squared is computed.
#' @return A data frame containing all measures of R-squared
#' @references
#' * `r2_d`: Royston P, Sauerbrei W. A new measure of prognostic separation in
#' survival data. Stat Med. 2004 Mar 15;23(5):723-48. doi: 10.1002/sim.1621.
#' PMID: 14981672.
#' * `r2_cs`, `rho_k`, `r2_pm`, `r2_r`, `rho_wa`: Royston, P. (2006). Explained
#' Variation for Survival Models. The Stata Journal, 6(1), 83–96.
#' https://doi.org/10.1177/1536867X0600600105
#' * `r2_xoq`: Ronghui Xu & John O' quigley (1999) A. R.2 type measure of
#' dependence for proportional hazards models, Journal of Nonparametric
#' Statistics, 12:1, 83-107, DOI: 10.1080/10485259908832799
#' @export
compute_rsquared <- function(data, f, cs_only = FALSE) {
  res <- data.frame(matrix(NA, nrow = 1, ncol = 0))

  model <- coxph(f, data)
  x2 <- 2 * (model$loglik[2] - model$loglik[1])
  res$r2_cs <- 1 - exp(-x2 / nrow(data))
  if (cs_only) {
    return(res)
  }
  sig_eps2 <- pi^2 / 6 # variance of the error term in Cox PH model
  rho_k <- 1 - exp(-x2 / sum(data$event))
  res$r2_r <- rho_k / (rho_k + sig_eps2 * (1 - rho_k))
  res$rho_k <- rho_k

  data <- data %>%
    mutate(
      cox_linear = model$linear.predictors,
      # -0.5 stabilizes (continuity correction)
      rankit = qnorm((rank(cox_linear) - 0.5) / n())
    )

  var_beta_x <- var(data$cox_linear)
  res$r2_pm <- var_beta_x / (var_beta_x + sig_eps2)
  res$rho_wa <- var_beta_x / (var_beta_x + 1)

  sigma2 <- (coxph(update.formula(f, ~ rankit), data)$coefficients[[1]])^2
  res$r2_d <- sigma2 / (sigma2 + sig_eps2)
  # TODO: add ref for r2_i?
  res$r2_i <- sigma2 / (sigma2 + 1)

  # xu1999 eq. (3.5)
  km <- survfit(Surv(time, event) ~ 1, data = data)
  # keep only time points with events and jump[i] = km$surv[i - 1] - km$surv[i]
  jump <- abs(diff(c(1, km$surv[which(km$n.event != 0)])))
  inter_val <- data %>%
    # `aeqSurv` gives the adjusted time points used by `survfit`
    mutate(time_adjusted = aeqSurv(Surv(time, event))[, 1]) %>%
    arrange(time) %>%
    # `rev(cumsum(rev(A)))` calculates `sum(df$A[i:nrow(df)])`
    mutate(
      # Effectively gives the number of patients as risk after combined with
      # `distinct` later
      Ci = n() + 1 - row_number(),
      Bi = rev(cumsum(rev(exp(cox_linear)))),
      # TODO: verify the correctness in multi-covariate case
      Ai = rev(cumsum(rev(cox_linear * exp(cox_linear)))),
    ) %>%
    filter(event != 0) %>%
    distinct(time_adjusted, .keep_all = TRUE)

  gamma <- 2 * sum(jump * (
    inter_val$Ai / inter_val$Bi - log(inter_val$Bi / inter_val$Ci)
  ))
  # Divided by `sum(jump)` to account for finite time, xu1999 eq. (3.6)
  res$rho_xoq <- 1 - exp(-gamma / sum(jump))
  return(res)
}
