library(survival)
library(dplyr)
library(furrr)
library(future)
library(purrr)
library(tidyr)

library(CovadjustSim)

set.seed(42)
strategy <- multisession
if (parallelly::supportsMulticore()) {
  strategy <- multicore
}
n_cores <- parallelly::availableCores()
initial_plan <- plan(strategy, workers = n_cores)
# Step 1.
# Set data frame for parameters, each row represents a simulation configuration.
params_grid <- expand_grid(
  shape = c(0.5, 1.0, 1.5),
  dropout = c(0.01, 0.1),
  t_cutoff = 5,
  hr_trt = c(0.4, 0.7),
  cindex = round(seq(0.55, 0.85, 0.1), 2),
  lambda = round(seq(0.1, 0.9, 0.1), 1)
) %>%
  # Set parameters `beta` and `kappa` for `simu_ph_weibull`
  fix_dist_params(
    beta_list = c(0.001, seq(0.1, 3, 0.1)),
    kappa_list = seq(-7.6, 4, 0.2),
    n_samples = 10000
  ) %>%
  arrange(t_cutoff, shape, desc(hr_trt), dropout, cindex, lambda) %>%
  rename(
    cindex0 = cindex,
    lambda0 = lambda
  ) %>%
  mutate(id = row_number())

# Two formulae / groups of covariates to be studied in each simulation
f1 <- Surv(time, event) ~ 1
f2 <- Surv(time, event) ~ X

# Step 2.
# For each simulation and for each formula, do a quick skim through the number
# of samples in the simulation and get the corresponding statistical power of
# the treatment effect estimation. Interpolate to get the range of sample sizes
# that cover the 80% threshold of interest.
# parallel across `params_grid`
df_range <- future_map_dfr(seq_len(nrow(params_grid)),
  function(i) {
    p <- params_grid %>% filter(id == i)
    simu <- simu_ph_weibull(
      p$shape, p$hr_trt, p$dropout, p$t_cutoff, p$beta, p$kappa
    )
    # At the beginning of the skim with small sample sizes, `survival::coxph`
    # may get some warnings on the convergence, which is of no concern as we are
    # only doing a rough estimation here.
    df1 <- power_skim(simu, f1) %>% suppressWarnings()
    df2 <- power_skim(simu, f2) %>% suppressWarnings()
    res <- c()
    res[c("n50_1", "n90_1")] <- interpolate_power(df1, c(0.5, 0.9))$n_samples
    res[c("n50_2", "n90_2")] <- interpolate_power(df2, c(0.5, 0.9))$n_samples
    return(res)
  },
  .options = furrr_options(seed = TRUE, chunk_size = 1)
)
params_grid <- params_grid %>% bind_cols(df_range)

# Step 3.
# For each simulation and each formula, estimate the sample size required to
# reach 80% statistical power by interpolation with the `power ~ n_samples`
# curve. Each curve is fitted by `n_points` points, at each point, `n_reps`
# repetitions are done to estimate the statistical power.
n_points <- 6
n_reps <- 10000
df_r2_obs_80 <- map_dfr(seq_len(nrow(params_grid)),
  function(i) {
    p <- params_grid %>% filter(id == i)
    simu <- simu_ph_weibull(
      p$shape, p$hr_trt, p$dropout, p$t_cutoff, p$beta, p$kappa
    )
    # parallel across `n_reps`
    res <- simu %>% estimate_sample_size_reduction(
      formula1 = f1, begin1 = p$n50_1, end1 = p$n90_1,
      formula2 = f2, begin2 = p$n50_2, end2 = p$n90_2,
      n_reps = n_reps, n_points = n_points
    ) %>%
      # cindex1 is always equal to 0.5, remove
      select(!cindex1) %>%
      rename(cindex = cindex2) %>%
    return(res)
  }
)

# Step 4.
# For each simulation, estimate various measures of r-squared for the formula
# `f`. Each estimation is the average over `n_reps_r2` repetitions and each
# repetition is based on a generated dataset of only the control arm and of
# `n_samples_r2` samples.
f <- Surv(time, event) ~ X
n_samples_r2 <- 1000
n_reps_r2 <- 1000
# parallel across `params_grid`
df_r2 <- future_map_dfr(seq_len(nrow(params_grid)),
  function(i) {
    p <- params_grid %>% filter(id == i)
    simu <- simu_ph_weibull(
      # `hr_trt = 1` is equivalent to only control arm
      p$shape, hr_trt = 1, p$dropout, p$t_cutoff, p$beta, p$kappa
    )
    res <- map_dfr(1:n_reps_r2, function(j) {
      return(simu %>%
        generate_trial_data(n_samples_r2) %>%
        compute_rsquared(f))
    }) %>%
      summarise_all(mean)
    return(res)
  },
  .options = furrr_options(seed = TRUE)
)

plan(initial_plan)  # restore future setting

parametric_study <- params_grid %>%
  bind_cols(df_r2_obs_80, df_r2) %>%
  relocate(
    id, shape, dropout, t_cutoff, hr_trt, cindex0, cindex, lambda0, lambda,
    beta, kappa, r2_obs_80
  )

usethis::use_data(parametric_study, overwrite = TRUE)
