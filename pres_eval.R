# Load required libraries
library(tidyverse)

## ============================================================================

library(EpiEstim)
library(EpiNow2)

library(incidence)
library(outbreak)

incidence_df <- EpiNow2::example_confirmed

options(mc.cores = 4)

## Just Looking at the historical period, and ignoring right-truction
## because EpiEstim doesn't do now-casting and thus

##
fixed_gamma <- Gamma(mean = 3, sd = 1, max = 10)

## SERIAL INTERVAL
generation_time <- Gamma(
  shape = Normal(9, 2.5), rate = Normal(3, 1.4), max = 10
)

## INCUBATION AND REPORTING DELAY
incubation_period <- LogNormal(
  meanlog = Normal(1.6, 0.05),
  sdlog = Normal(0.5, 0.05),
  max = 14
)

reporting_delay <- LogNormal(meanlog = 0.5, sdlog = 0.5, max = 10)

delay <- incubation_period + reporting_delay

##
rt_prior <- LogNormal(mean = 2, sd = 1)


##
estRt <- epinow(
  example_confirmed,
  generation_time = gt_opts(generation_time),
  delays = delay_opts(delay),
  rt = rt_opts(prior = rt_prior)
)
