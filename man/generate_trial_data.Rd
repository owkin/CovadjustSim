% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_trial_data.R
\name{generate_trial_data}
\alias{generate_trial_data}
\alias{generate_trial_data.simu_ph_weibull}
\title{Generate a dataset simulating a clinical trial}
\usage{
generate_trial_data(x, ...)

\method{generate_trial_data}{simu_ph_weibull}(x, n_samples, ...)
}
\arguments{
\item{x}{An object of a subclass of \link{simu_ph_base}.}

\item{...}{further arguments passed to other methods.}

\item{n_samples}{Number of samples in the returned data frame.}
}
\value{
A data frame of simulated clinical trial with at least three columns:
time: numeric, the survival time.
event: integer 1 or 0, indicating either an event (e.g. death) or right
censoring.
treatment: integer 1 or 0, indicating whether the patient is in the
treatment arm or in the control arm.
Additionally, other columns can exist which are considered as potential
covariates in the study.
}
\description{
No default method. Each class (subclass of \link{simu_ph_base}) should implement
this method.
}
