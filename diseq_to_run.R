library(data.table)
library(magrittr)
library(parallel)
library(ggplot2)

#################
# Preliminaries #
#################

setwd('X:/PSF/_Common/Személyek/BanaiÁdám/Horváth Ákos/panel/diseq_R')

source('02_Kodok/diseq_maxim.R')

dt <- fread('01_Adatbazisok/01_panel_diseq.csv', integer64 = 'double')

# Sample restriction to relevant firms
dt[, ind_factor := factor(ind)]
dt_model <- dt[d12_apply == 1 | d13_apply == 1 | ydis13_newloans > 0]

set.seed(1984)
dt_model[, r := runif(.N)]

########################
# Estimating the model #
########################

fixed_noncons <- c(0, 1, 1, 1, 1,
                   0, 1, 1, 1, 1, 1,
                   0, 0)
fixed_cons <- c(1, 0, 0, 0, 0
                1, 0, 0, 0, 0, 0,
                0, 0)

diseq_model_cens_nocorr_simple_iter <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha,
    supply_formula = ydis13_newloans ~ x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears,
    data = dt_model[r <= 0.05], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 5000),
    lb = c(-1000, -100, -100, -100, -100,
           -1000, -100, -100, -100, -100, -100,
           -5, -5),
    ub = c(1000, 100, 100, 100, 100,
           1000, 100, 100, 100, 100, 100,
           5, 5),
    init = c(0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0,
             0, 0)
    fixed_params = fixed_noncons
  ) %>%
  refitdiseq(
    fixed_params = fixed_cons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = fixed_noncons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = fixed_cons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = fixed_noncons,
    continue = TRUE,
    initpop = FALSE
  ) %>%
  refitdiseq(
    fixed_params = fixed_cons,
    continue = TRUE,
    initpop = FALSE
  )

diseq_model_cens_nocorr_simple_iter %>%
  export_all(folder = '03_Eredmenyek/simple_DE_2013_5pc_iter')

###############
# Diagnostics #
###############

print(summary(diseq_model_cens_nocorr_simple_iter))

optim_trace_cons <- diseq_model_cens_nocorr_bovebb_de$optim_trace[c(1, 3, 5)]
optim_trace_noncons <- diseq_model_cens_nocorr_bovebb_de$optim_trace[c(1, 3, 5)]

loglik_dt_ind <- diseq_model_cens_nocorr_bovebb_de$optim_trace %>%
  lapply(function (trace) as.data.table(trace[[1]])) %>%
  rbindlist()
loglik_dt_ind[, iter_no := .I]
ggplot(loglik_dt_ind) + geom_line(aes(x = iter_no, y = V1))

params_dt_ind_cons <- gyumi_model_DE_iter$optim_trace[c(1, 3, 5)] %>%
  lapply(function (trace) as.data.table(trace[[2]])) %>%
  rbindlist() %>%
  .[, iter_no := .I] %>%
  melt(id.vars = "iter_no")
ggplot(params_dt_ind_cons) + geom_line(aes(x = iter_no, y = value, color = variable))

params_dt_ind_noncons <- gyumi_model_DE_iter$optim_trace[c(2, 4, 6)] %>%
  lapply(function (trace) as.data.table(trace[[2]])) %>%
  rbindlist() %>%
  .[, iter_no := .I] %>%
  melt(id.vars = "iter_no")
ggplot(params_dt_ind_noncons) + geom_line(aes(x = iter_no, y = value, color = variable))
