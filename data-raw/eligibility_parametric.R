library(dplyr)
library(purrr)
library(survival)
library(withr)

library(CovadjustSim)

set.seed(1)
# Define shared simulation parameters for each simulation
params_common <- data.frame(
  shape = 0.5,
  hr_trt = 0.7,
  dropout = 0,
  t_cutoff = 5,
  cindex = 0.65,
  lambda = 0.9
) %>%
  # Set values of `beta` and `kappa` that give the desired `cindex` and `lambda`
  fix_dist_params(n_samples = 10000, n_cores = 8) %>%
  select(!c(cindex, lambda))

# Two sets of covariates
formulae <- list(
  Unadjusted = Surv(time, event) ~ 1,
  Adjusted = Surv(time, event) ~ X
)
# Two eligibility criteria
params1 <- params_common %>% mutate(inclusion = 1.0)
simu1 <- do.call(simu_ph_weibull, params1)
params2 <- params_common %>% mutate(inclusion = 0.8)
simu2 <- do.call(simu_ph_weibull, params2)

# Fix (roughly) the number of events in each simulation so that in the final
# power ~ n_events curve, different eligibility criteria have aligned points.
# Here the less restrictive and unadjusted model is taken as reference.
df_power <- power_skim(simu1, formulae$Unadjusted, n_reps = 1000)
# The range of the number of events corresponds to 50% and 90% statistical power
n_events_range <- interpolate_power(df_power, c(0.5, 0.9))$n_events
n_events_list <- seq(n_events_range[1], n_events_range[2], length.out = 9)

# Wrapper function for `compute_power()`
func <- function(n_samples, simu, n_reps, formulae) {
  return(compute_power(simu, n_reps, n_samples, formulae,
    get_cindex = TRUE, get_lambda = TRUE, n_cores = parallelly::availableCores()
  ))
}

n_reps <- 10000
# Run simulaitons
event_freq1 <- mean(sapply(1:1000, function(i) {
  d <- simu1 %>% generate_trial_data(1000)
  return(sum(d$event) / nrow(d))
}))
n_samples_list1 <- ceiling(n_events_list / event_freq1)
res1 <- map_dfr(n_samples_list1, func, simu1, n_reps, formulae)

event_freq2 <- mean(sapply(1:1000, function(i) {
  d <- simu2 %>% generate_trial_data(1000)
  return(sum(d$event) / nrow(d))
}))
n_samples_list2 <- ceiling(n_events_list / event_freq2)
res2 <- map_dfr(n_samples_list2, func, simu2, n_reps, formulae)

# Gather results
eligibility_parametric <- bind_rows(
  res1 %>% mutate(inclusion = "100%"),
  res2 %>% mutate(inclusion = "80%")
) %>%
  mutate(inclusion = factor(inclusion, levels = c("100%", "80%")))

usethis::use_data(eligibility_parametric, overwrite = TRUE)
