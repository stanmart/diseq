library(data.table)
library(magrittr)
library(parallel)
library(ggplot2)

#################
# Preliminaries #
#################

source('diseq_maxim.R')

dt <- fread('01_Adatbazisok/01_panel_diseq.csv', integer64 = 'double')

# Sample restriction to relevant firms
dt[, ind_factor := factor(ind)]
dt_model <- dt[d12_apply == 1 | d13_apply == 1 | ydis13_newloans > 0]

set.seed(1984)
dt_model[, r := runif(.N)]

########################
# Estimating the model #
########################

model_ols <- lm(
  formula = ydis13_newloans ~ x12_age + x12_size + x12_trade + x12_soa + x12_wcratio +
    x12_liqasset + x12_growtha + x12_growths + c12_own,
  data = dt_model[r <= 0.05 & ydis13_newloans > 0]
)

model_iwls <- glm(
  formula = ydis13_newloans ~ x12_age + x12_size + x12_trade + x12_soa + x12_wcratio +
    x12_liqasset + x12_growtha + x12_growths + c12_own,
  family = gaussian(link = "identity"),
  data = dt_model[r <= 0.05 & ydis13_newloans > 0]
)

model_diseq_local <- fitdiseq(
  demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_trade + x12_soa + x12_wcratio +
    x12_liqasset + x12_growtha + x12_growths + c12_own,
  supply_formula = ydis13_newloans ~ 1,
  data = dt_model[r <= 0.05 & ydis13_newloans > 0],
  optimizer = 'optim',
  lb = c(-1000, -100, -100, -100, -100, -100, -100, -100, -100, -100,
         -Inf,
         -5, -Inf),
  ub = c(1000, 100, 100, 100, 100, 100, 100, 100, 100, 100,
         Inf,
         5, Inf),
  init = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           Inf,
           1, 1),
  fixed_params = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   1,
                   0, 1),
  likelihood_bias = 0
)

model_diseq_global <- fitdiseq(
  demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_trade + x12_soa + x12_wcratio +
    x12_liqasset + x12_growtha + x12_growths + c12_own,
  supply_formula = ydis13_newloans ~ 1,
  data = dt_model[r <= 0.05 & ydis13_newloans > 0],
  optimizer = 'DE',
  control = DEoptim.control(trace = 100, itermax = 5000),
  lb = c(-1000, -100, -100, -100, -100, -100, -100, -100, -100, -100,
         -Inf,
         -5, -Inf),
  ub = c(1000, 100, 100, 100, 100, 100, 100, 100, 100, 100,
         Inf,
         5, Inf),
  init = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           Inf,
           1, 1),
  fixed_params = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   1,
                   0, 1),
  likelihood_bias = 0
)

model_diseq_local %>% export_all(folder = '03_Eredmenyek/supply_switched_off_local_0604')
model_diseq_global %>% export_all(folder = '03_Eredmenyek/supply_switched_off_global_0604')

estimates <- rbind(
  model_ols$coefficients,
  model_iwls$coefficients,
  model_diseq_local$coefficients[model_diseq_local$coef_indices$beta_demand],
  model_diseq_global$coefficients[model_diseq_global$coef_indices$beta_demand]
)

reg_table <- estimates %>% 
  t() %>% 
  as.data.table() %>%
  setnames(c("OLS", "IWLS", "diseq local", "diseq_global")) %>% 
  .[, variable := names(model_ols$coefficients)] %>% 
  setcolorder("variable") %>% 
  knitr::kable() %>% 
  print()

