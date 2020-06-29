source('diseq_maxim.R')
library(data.table)

N <- 2000

beta1 <- 1
beta2 <- -0.3
beta3 <- -0.7
beta4 <- 0.8
rho <- 0.5

dt <- data.table(
  alma = rnorm(N, 1, 1),
  korte = rnorm(N, 2, 1),
  szilva = rnorm(N, 3, 1),
  e1 = rnorm(N)
)
dt[, e2 := rho * e1 + sqrt(1 - rho^2) * rnorm(N)]

dt[, gyumi_kereslet := beta1 * alma + beta2 * szilva + e1]
dt[, gyumi_kinalat := beta3 * korte + beta4 * szilva + e2]

dt[, gyumi := pmin(gyumi_kereslet, gyumi_kinalat)]
dt[, gyumi_poz := pmax(0, gyumi)]

gyumi_model_SA <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + szilva,
  supply_formula = gyumi_poz ~ korte + szilva,
  data = dt,
  lb = rep(-5, 9),
  ub = rep(5, 9),
  optimizer = 'SA',
  corr = TRUE
)

dt[, truncated := gyumi_poz == 0]
dt[, prob_uncond := predict(gyumi_model_SA, type = "prob_supply_constrained", conditional = FALSE)]
dt[, prob_cond := predict(gyumi_model_SA, type = "prob_supply_constrained", conditional = TRUE)]
dt[, exp_deficit := predict(gyumi_model_SA, type = "expected_deficit", conditional = TRUE)]

ggplot(dt) + geom_point(aes(x = prob_uncond, y = prob_cond, color = truncated))
ggplot(dt) + geom_point(aes(x = prob_uncond, y = exp_deficit, color = truncated))
ggplot(dt) + geom_point(aes(x = prob_cond, y = exp_deficit, color = truncated))
