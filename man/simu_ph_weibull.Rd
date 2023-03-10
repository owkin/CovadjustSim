% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{simu_ph_weibull}
\alias{simu_ph_weibull}
\title{S3 class of trial simulation based on Weibull and proportional hazard model}
\usage{
simu_ph_weibull(shape, hr_trt, dropout, t_cutoff, beta, kappa, inclusion = 1)
}
\arguments{
\item{shape}{The shape parameter of the Weibull distribution.}

\item{hr_trt}{The treatment hazard ratio.}

\item{dropout}{The rate parameter in the exponential distribution to model
the right censoring.}

\item{t_cutoff}{The cut-off time (in year) beyond which all patients will be
right censored.}

\item{beta, kappa}{Parameters of the hazard term \code{kappa + beta * X} where
\code{X} is the covariate. Their values are often fixed so as to reach a desired
cumulative incidence \code{lambda} and a desired concordance \verb{cindex' of the fitted proportional hazard model at }t_cutoff` in the control arm.}

\item{inclusion}{A numerical value between ]0, 1] to simulate the inclusion
criterion. Specifically, it is a quantile threshold of the covariate \code{X},
only samples with a covariate value below the threshold are included in the
trial.}
}
\value{
An object of class \code{c(simu_ph_weibull, simu_ph_base)}
}
\description{
Baseline survival times are modeled by the Weibull distribution.
Covariate is a variable \code{X} following the standard normal distribution.
Overall survival is modeled by the proportional hazard model.
Right censoring time is modeled by an exponential distribution with specified
dropout rate.
}
\seealso{
\link{simu_ph_base}
}
