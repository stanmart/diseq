library(data.table)
library(magrittr)
library(parallel)
library(ggplot2)

#################
# Preliminaries #
#################

setwd('X:/PSF/_Common/Személyek/BanaiÁdám/Horváth Ákos/panel/diseq_R')

source('02_Kodok/diseq-to_run_v1/diseq_maxim.R')

dt <- fread('01_Adatbazisok/01_panel_diseq.csv', integer64 = 'double')

# Sample restriction to relevant firms
dt[, ind_factor := factor(ind)]
dt_model <- dt[d12_apply == 1 | d13_apply == 1 | ydis13_newloans > 0]

set.seed(1984)
dt_model[, r := runif(.N)]

########################
# Estimating the model #
########################

optim_cons <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1)
optim_noncons <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   1, 1)
optim_stdev <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 0, 0)

diseq_model_cens_nocorr_simple_iter <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_trade + x12_soa + x12_wcratio +
      x12_liqasset + x12_growtha + x12_growths + c12_own,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince,
    data = dt_model[r <= 0.05],
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 5000),
    lb = c(-1000, -100, -100, -100, -100, -100, -100, -100, -100, -100,
           -1000, -100, -100, -100, -100, -100, -100, -100, -100, -100,
           -5, -5),
    ub = c(1000, 100, 100, 100, 100, 100, 100, 100, 100, 100,
           1000, 100, 100, 100, 100, 100, 100, 100, 100, 100,
           5, 5),
    init = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             1, 1),
    fixed_params = optim_noncons,
    likelihood_bias = exp(-10)
  ) %>%
  refitdiseq(
    fixed_params = optim_cons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = optim_stdev,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = optim_noncons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = optim_cons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = optim_stdev,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = optim_noncons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = optim_cons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = optim_stdev,
    continue = TRUE,
    initpop = FALSE
  )

diseq_model_cens_nocorr_simple_iter %>%
  export_all(folder = '03_Eredmenyek/simple_DE_2013_5pc_iter_akos_0416')

###############
# Diagnostics #
###############

print(summary(diseq_model_cens_nocorr_simple_iter))

optim_trace_noncons <- diseq_model_cens_nocorr_simple_iter$optim_trace[c(1, 4, 7)]
optim_trace_cons <- diseq_model_cens_nocorr_simple_iter$optim_trace[c(2, 5, 8)]
optim_trace_stdev <- diseq_model_cens_nocorr_simple_iter$optim_trace[c(3, 6, 9)]

coefnames <- names(diseq_model_cens_nocorr_simple_iter$coefficients)

coefnames_cons <- coefnames[!optim_cons]
coefnames_noncons <- coefnames[!optim_noncons]
coefnames_stdev <- coefnames[!optim_stdev]

coefnames_noncons_demand <- intersect(
  coefnames_noncons, 
  coefnames[diseq_model_cens_nocorr_simple_iter$coef_indices$beta_demand]
)
coefnames_noncons_supply <- intersect(
  coefnames_noncons, 
  coefnames[diseq_model_cens_nocorr_simple_iter$coef_indices$beta_supply]
)

loglik_dt_ind <- diseq_model_cens_nocorr_simple_iter$optim_trace %>%
  lapply(function (trace) as.data.table(trace[[1]])) %>%
  rbindlist()
loglik_dt_ind[, iter_no := .I]
ggplot(loglik_dt_ind) + geom_line(aes(x = iter_no, y = V1))

params_dt_ind_cons <- optim_trace_cons %>%
  lapply(function (trace) as.data.table(trace[[2]])) %>%
  rbindlist() %>%
  setnames(coefnames_cons) %>%
  .[, iter_no := .I] %>%
  melt(id.vars = "iter_no")
ggplot(params_dt_ind_cons) +
  geom_line(aes(x = iter_no, y = value, color = variable)) +
  ggtitle("Constants")

params_dt_ind_noncons <- optim_trace_noncons %>%
  lapply(function (trace) as.data.table(trace[[2]])) %>%
  rbindlist() %>%
  setnames(coefnames_noncons) %>%
  .[, iter_no := .I] %>%
  melt(id.vars = "iter_no")
ggplot(params_dt_ind_noncons[variable %in% coefnames_noncons_demand]) +
  geom_line(aes(x = iter_no, y = value, color = variable)) +
  ggtitle("Non-constamts (demand side)")
ggplot(params_dt_ind_noncons[variable %in% coefnames_noncons_supply]) +
  geom_line(aes(x = iter_no, y = value, color = variable)) +
  ggtitle("Non-constamts (supply side)")

params_dt_ind_stdev <- optim_trace_stdev %>%
  lapply(function (trace) as.data.table(trace[[2]])) %>%
  rbindlist() %>%
  setnames(coefnames_stdev) %>%
  .[, iter_no := .I] %>%
  melt(id.vars = "iter_no")
ggplot(params_dt_ind_noncons) +
  geom_line(aes(x = iter_no, y = value, color = variable)) +
  ggtitle("Dispersion parameters")
