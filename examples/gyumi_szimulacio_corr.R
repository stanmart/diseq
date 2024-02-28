source('diseq_maxim.R')
library(data.table)

N <- 10000

beta1 <- 1
beta2 <- -0.3
beta3 <- -0.7
beta4 <- 0.8
rho = 0.5

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

# Simulated annealing

t0 <- Sys.time()
gyumi_model_SA_nocorr <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + szilva,
  supply_formula = gyumi_poz ~ korte + szilva,
  data = dt,
  optimizer = 'SA',
  corr = FALSE
)
running_time_SA_nocorr <- difftime(Sys.time(), t0, units = 'min')
print(gyumi_model_SA_nocorr)

gyumi_summary_SA_nocorr <- summary(gyumi_model_SA_nocorr)
print(gyumi_summary_SA_nocorr)

t0 <- Sys.time()
gyumi_model_SA <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + szilva,
  supply_formula = gyumi_poz ~ korte + szilva,
  data = dt,
  optimizer = 'SA',
  corr = TRUE
)
running_time_SA <- difftime(Sys.time(), t0, units = 'min')
print(gyumi_model_SA)

gyumi_summary_SA <- summary(gyumi_model_SA)
print(gyumi_summary_SA)

# Differential evolution
t0 <- Sys.time()
gyumi_model_DE <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + szilva,
  supply_formula = gyumi_poz ~ korte + szilva,
  data = dt,
  optimizer = 'DE'
)
running_time_DE <- difftime(Sys.time(), t0, units = 'min')
print(gyumi_model_DE)

gyumi_summary_DE <- summary(gyumi_model_DE)
print(gyumi_summary_DE)

model_comparison <- rbind(cbind(gyumi_model_SA$coefficients, gyumi_model_DE$coefficients),
                          cbind(gyumi_model_SA$log_likelihood, gyumi_model_DE$log_likelihood),
                          c(running_time_SA, running_time_DE))
colnames(model_comparison) <- c('SA', 'DE')
rownames(model_comparison) <- c(rownames(model_comparison)[1 : (nrow(model_comparison) - 2)],
                                'Log-likelihood', 'Running time (min)')
model_comparison <- cbind(model_comparison, (model_comparison[, 'SA'] - model_comparison[, 'DE']) / model_comparison[, 'DE'] * 100)
colnames(model_comparison) <- c('SA', 'DE', 'Percentage diff.')
print(model_comparison)
