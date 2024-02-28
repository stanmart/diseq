source('diseq_maxim.R')
library(data.table)

N <- 10000

beta1 <- 1
beta2 <- -0.3
beta3 <- -0.7
beta4 <- 0.8

dt <- data.table(
  alma = rnorm(N, 1, 1),
  korte = rnorm(N, 2, 1),
  szilva = rnorm(N, 3, 1),
  e1 = rnorm(N, 0, 2),
  e2 = rnorm(N, 0, 2)
)

dt[, gyumi_kereslet := beta1 * alma + beta2 * szilva + e1]
dt[, gyumi_kinalat := beta3 * korte + beta4 * szilva + e2]

dt[, gyumi := pmin(gyumi_kereslet, gyumi_kinalat)]
dt[, gyumi_poz := pmax(0, gyumi)]

model_ols <- lm(
  formula = gyumi_poz ~ alma + korte,
  data = dt[gyumi_poz > 0]
)

model_iwls <- glm(
  formula = gyumi_poz ~ alma + korte,
  family = gaussian(link = "identity"),
  data = dt[gyumi_poz > 0]
)

model_diseq_local <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + korte,
  supply_formula = gyumi_poz ~ 1,
  init = c(0, 0, 0, Inf, 0, 1),
  fixed_params = c(0, 0, 0, 1, 0, 1),
  data = dt[gyumi_poz > 0],
  optimizer = "optim"
)

model_diseq_global <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + korte,
  supply_formula = gyumi_poz ~ 1,
  init = c(0, 0, 0, Inf, 0, 1),
  fixed_params = c(0, 0, 0, 1, 0, 1),
  data = dt[gyumi_poz > 0],
  optimizer = "SA"
)
