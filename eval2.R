
library(EpiEstim)
library(EpiNow2)
library(ggplot2)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# AGGREGATED INCIDENCE DATA
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

incidence_df <- EpiNow2::example_confirmed
incidence_df <- incidence_df[10:50, ]

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# GENERATION TIME and SERIAL INTERVAL
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

generation_interval_pmf <- c(
  0.0610, 0.1540, 0.2198, 0.2178, 0.1536, 0.1122, 0.0486, 0.0224,
  0.0078, 0.0022, 0.0004, 0.0002
)

## Has to start with 0 !
gi_pmf <- NonParametric(pmf = c(0, generation_interval_pmf))

## ASSUME THAT SERIAL INTERVAL = GENERATION INTERVAL
## also has to start with 0
serial_interval_pmf <- c(0, generation_interval_pmf)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# DELAYS
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# 1) Infection to taking a test
infect_to_test_pmf <- c(
  0, 0.0001, 0.0030, 0.1422, 0.2714, 0.2664,
  0.1832, 0.0846, 0.0340, 0.0136, 0.0014, 0.0001
)

infect_to_test <- NonParametric(pmf = infect_to_test_pmf)


# 2) time from taking a test to it getting into a state database
reporting_delay_pmf <- c(
  0.3786, 0.3724, 0.1662, 0.0622, 0.0166, 0.0034, 0.0006
)

sym_report_delay_pmf <- NonParametric(pmf = reporting_delay_pmf)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# ADD RIGHT TRUNCATION TO DATA
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

tail_Px <- rev(cumsum(reporting_delay_pmf))
px_len <- length(tail_Px)

incidence_df$Px = 1
cases_len <- nrow(incidence_df)

incidence_df$Px[(cases_len - px_len + 1):cases_len] <- tail_Px

incidence_df$confirm_true <- incidence_df$confirm

incidence_df$confirm <- incidence_df$confirm * incidence_df$Px

incidence_df$confirm <- as.integer(incidence_df$confirm)

plot(x = incidence_df$date,
     y = incidence_df$confirm,
     type = 'l')

lines(x = incidence_df$date,
      y = incidence_df$confirm_true,
      type = 'l', col = 'red')


# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# ----------------------------------------------------------------------------

get_Rt <- function(m, w) {

  R <- rep(NA, length(m))

  for(t in 2:length(m)) {
    tau_end = min(length(w), t - 1)
    c_mat <- w[1:tau_end] %*% m[t - 1:tau_end]
    R[t] <- m[t] / c_mat
  }

  return(R)
}

#
Rt_df <- data.frame(report_date = incidence_df$date)
Rt_df$model <- 'Rt_calc'

Rt_df$Rt <- get_Rt(incidence_df$confirm, generation_interval_pmf)

#
INFECT_TO_TEST_SHIFT = round(weighted.mean(x = 1:length(infect_to_test_pmf),
                                       w = infect_to_test_pmf))

REPORTINGDELAY_SHIFT = round(weighted.mean(x = 1:length(reporting_delay_pmf),
                                           w = reporting_delay_pmf))

Rt_df$date = Rt_df$report_date -  INFECT_TO_TEST_SHIFT - REPORTINGDELAY_SHIFT
Rt_df$date_int <- 1:nrow(Rt_df) - INFECT_TO_TEST_SHIFT - REPORTINGDELAY_SHIFT

Rt_df$report_date <- NULL

head(Rt_df)

Rt_df$Rt_lb <- NA
Rt_df$Rt_ub <- NA

ggplot(Rt_df) +
  geom_line(aes(x = date_int, y = Rt))

head(Rt_df)
# ----------------------------------------------------------------------------

##
epiestim_df <- incidence_df

epiestim_df$I <- epiestim_df$confirm

# Estimate R DAILY
getR <- EpiEstim::estimate_R(
  incid = epiestim_df ,
  method = "non_parametric_si",
  config = make_config(list(
    si_distr = serial_interval_pmf
  ))
)

## assign where between t_start and t_end is dt_int
## has to be an integer
## OBVIOUSLY THIS CAN HAVE A HUGE IMPACT ON YOUR ESTIMATES
getR$R$date_int <- round((getR$R$t_end - getR$R$t_start) / 2 + getR$R$t_start)
getR$R$date_int <- round((getR$R$t_end - getR$R$t_start) / 2 + getR$R$t_start)

EpiEstim_R <- getR$R[, c('date_int', 'Median(R)',
                       'Quantile.0.05(R)', 'Quantile.0.95(R)')]

names(EpiEstim_R) <- c("date_int", "Rt", "Rt_lb", "Rt_ub")

EpiEstim_R$date <- NA
EpiEstim_R$model <- 'EpiEstim'

## shift the dates
EpiEstim_R$date_int = EpiEstim_R$date_int -
  INFECT_TO_TEST_SHIFT -
  REPORTINGDELAY_SHIFT

##
head(EpiEstim_R)

Rt_df <- rbind(Rt_df, EpiEstim_R)


ggplot(Rt_df) +
  geom_pointline(aes(x = date_int, y = Rt, color = model))

# ----------------------------------------------------------------------------
res_epinow <- epinow(
  data            = incidence_df[, c('date', 'confirm')],
  generation_time = generation_time_opts(gi_pmf),
  delays          = delay_opts(infect_to_test),
  truncation      = trunc_opts(sym_report_delay_pmf),
  backcalc        = backcalc_opts(prior = 'reports'),
  # rt = rt_opts(rw = 1),
  stan            = stan_opts(chains = 4, cores = 4)
)

R_df <- subset(res_epinow$estimates$summarised, variable == 'R')

R_df <-




#
INCUBATION_SHIFT = round(weighted.mean(x = all_data$incubation$Day,
                                       w = all_data$incubation$Px))

REPORTINGDELAY_SHIFT = round(weighted.mean(x = all_data$reporting_delay$Day,
                                           w = all_data$reporting_delay$Px))

