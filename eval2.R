
library(EpiEstim)
library(EpiNow2)
library(ggplot2)
library(ggpubr)
library(patchwork)

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
# ADD RIGHT TRUNCATION TO DATA BASED ON REPORTING DELAY
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

incidence_df$date_int <- 1:nrow(incidence_df) - 1

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
Rt_df$date_int <- 1:nrow(Rt_df) - 1 - INFECT_TO_TEST_SHIFT - REPORTINGDELAY_SHIFT

Rt_df$report_date <- NULL

head(Rt_df)

Rt_df$Rt_lb <- NA
Rt_df$Rt_ub <- NA

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
getR$R$date_int <- getR$R$t_end

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

Rt_df2 <- rbind(Rt_df, EpiEstim_R)

## shift the dates
EpiEstim_R$date_int = EpiEstim_R$date_int - mean((getR$R$t_end - getR$R$t_start) / 2)
EpiEstim_R$model <- 'EpiEstim (midpoint)'
##
head(EpiEstim_R)

Rt_df3 <- rbind(Rt_df2, EpiEstim_R)
head(Rt_df3)

# ----------------------------------------------------------------------------
res_epinow <- epinow(
  data            = incidence_df[, c('date', 'confirm')],
  generation_time = generation_time_opts(gi_pmf),
  delays          = delay_opts(infect_to_test),
  truncation      = trunc_opts(sym_report_delay_pmf),
  backcalc        = backcalc_opts(prior = 'reports'),
  stan            = stan_opts(chains = 4, cores = 4)
)

plot(res_epinow)

R_df <- subset(res_epinow$estimates$summarised, variable == 'R')
R_df$date_int <- 1:nrow(R_df) - 1
R_df$model <- 'EpiNow2'
R_df$Rt <- R_df$median
R_df$Rt_lb <- R_df$lower_90
R_df$Rt_ub <- R_df$upper_90

cases_df <- subset(res_epinow$estimates$summarised, variable == 'reported_cases')
cases_df$date_int <- 1:nrow(R_df) - 1
cases_df$model <- 'EpiNow2'
cases_df$reports <- cases_df$median
cases_df$reports_lb <- cases_df$lower_90
cases_df$reports_ub <- cases_df$upper_90


EpiNow2_df <- R_df[, c('model', 'Rt', 'date', 'date_int', 'Rt_lb', 'Rt_ub')]

Rt_df4 <- rbind(Rt_df3, EpiNow2_df)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 4
cols = gg_color_hue(n)
cols

col_EpiEstim <- cols[1]
col_EpiNow2 <- cols[3]

plot_rt1 <- ggplot(subset(Rt_df4, date_int >= 0 &
                            ! (model %in% c('Rt_calc', 'EpiEstim (midpoint)'))))+
  theme_classic2() +
  geom_hline(yintercept = 1, linetype = '11') +
  geom_vline(xintercept = c(40 + 1.5 - length(reporting_delay_pmf), 40.5 ),
             linetype = '41') +
  geom_ribbon(aes(x = date_int, ymin = Rt_lb, ymax = Rt_ub, fill = model),
              alpha = 0.25) +
  geom_line(aes(x = date_int, y = Rt, color = model)) +
  scale_color_manual(name = 'Model', values = c(col_EpiEstim,col_EpiNow2 )) +
  scale_fill_manual(name = 'Model', values = c(col_EpiEstim,col_EpiNow2 )) +
  coord_cartesian(xlim = c(0, 50)) + ylab("R(t)") +
  annotate('text', x = 10, y = 2.5, label = 'Historical period') +
  annotate('text', x = 37.5, y = 2.5, label = 'Nowcast') +
  annotate('text', x = 43.5, y = 2.5, label = 'Forecast') +
  xlab(NULL)

plot_report <- ggplot(incidence_df) + theme_classic2() +
  geom_col(aes(x = date_int, y = confirm), fill = grey(0.75), width = 0.5) +
  geom_vline(xintercept = c(40 + 1.5 - length(reporting_delay_pmf), 40.5 ),
             linetype = '41') +
  coord_cartesian(xlim = c(0, 50)) + ylab("Reported cases") +
  geom_ribbon(data = cases_df, aes(x = date_int,
                                   ymin = reports_lb, ymax = reports_ub),
              alpha = 0.25, fill = col_EpiNow2) +
  geom_line(data = cases_df, aes(x = date_int, y = reports), color = col_EpiNow2) +
  xlab('Date') +
  annotate('text', x = 10, y = 8000, label = 'Historical period') +
  annotate('text', x = 37.5, y = 8000, label = 'Nowcast') +
  annotate('text', x = 43.5, y = 8000, label = 'Forecast')


plot_rt1 / (plot_report) + plot_layout(guides = "collect")
dev.size()
ggsave('Eval.png', width = 6*1.5, height = 3.5*1.5)
