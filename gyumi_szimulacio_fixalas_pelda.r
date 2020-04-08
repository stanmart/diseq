source('diseq_maxim.R')
library(data.table)
library(magrittr)

N <- 5000

alpha1 <- 0.5
beta1 <- 1
beta2 <- -0.3
alpha2 <- 0
beta3 <- -0.7
beta4 <- 0.8

dt <- data.table(
  alma = rnorm(N, 1, 1),
  korte = rnorm(N, 2, 1),
  szilva = rnorm(N, 3, 1),
  e1 = rnorm(N, 0, 2),
  e2 = rnorm(N, 0, 2)
)

dt[, gyumi_kereslet := alpha1 + beta1 * alma + beta2 * szilva + e1]
dt[, gyumi_kinalat := alpha2 + beta3 * korte + beta4 * szilva + e2]

dt[, gyumi := pmin(gyumi_kereslet, gyumi_kinalat)]
dt[, gyumi_poz := pmax(0, gyumi)]

gyumi_model_SA <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + szilva,
  supply_formula = gyumi_poz ~ korte + szilva,
  data = dt,
  optimizer = 'SA',
  init = c(0, 0, 0, 0, 0, 0, 0, 0)
)

gyumi_model_SA_iter <- fitdiseq(
  demand_formula = gyumi_poz ~ alma + szilva,
  supply_formula = gyumi_poz ~ korte + szilva,
  data = dt,
  optimizer = 'SA',
  init = c(0, 0, 0, 0, 0, 0, 0, 0),
  fixed_params = c(1, 0, 0, 1, 0, 0, 0, 0)
) %>% refitdiseq(
  fixed_params = c(0, 1, 1, 0, 1, 1, 1, 1),
  continue = TRUE
) %>% refitdiseq(
  fixed_params = c(1, 0, 0, 1, 0, 0, 0, 0),
  continue = TRUE
) %>% refitdiseq(
  fixed_params = c(0, 1, 1, 0, 1, 1, 1, 1),
  continue = TRUE
)

model_comparison <- rbind(cbind(gyumi_model_SA$coefficients, gyumi_model_SA_iter$coefficients),
                          cbind(gyumi_model_SA$log_likelihood, gyumi_model_SA_iter$log_likelihood))
colnames(model_comparison) <- c('full', 'iterative')
rownames(model_comparison) <- c(rownames(model_comparison)[1 : (nrow(model_comparison) - 1)],
                                'Log-likelihood')
print(model_comparison)
