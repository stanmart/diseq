library(data.table)
library(magrittr)
library(parallel)
library(ggplot2)

setwd('X:/PSF/_Common/Személyek/BanaiÁdám/Horváth Ákos/panel/diseq_R')

source('02_Kodok/diseq_maxim.R')

dt <- fread('01_Adatbazisok/01_panel_diseq.csv', integer64 = 'double')

# Sample restriction to relevant firms
dt[, ind_factor := factor(ind)]

dt_model <- dt[d12_apply == 1 | d13_apply == 1 | ydis13_newloans > 0]

set.seed(1984)
dt_model[, r := runif(.N)]

# Model with no correlation, but with censoring - first without industry dummies

diseq_model_cens_nocorr <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.05], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 5000),
    lb = c(-1000, -10, -10, -10, -10, -10, -10, -10, -10,# rep(-5, 13),
           -1000, -10, -10, -10, -10, -10, -10, -10, -10, -10,# rep(-5, 13),
           -5, -5),
    ub = c(1000, 10, 10, 10, 10, 10, 10, 10, 10,# rep(5, 13),
           1000, 10, 10, 10, 10, 10, 10, 10, 10, 10,# rep(5, 13),
           5, 5)
  ) %>% export_all(folder = '03_Eredmenyek/masodik_proba_DE_2013_5pc') %>% 
  refitdiseq(
    optimizer = 'SA',
    control = list('verbose' = TRUE, 'max.time' = 15*60*60),
    lb = c(-1000, -100, -100, -100, -100, -100, -100, -100, -100,# rep(-5, 13),
           -1000, -100, -100, -100, -100, -100, -100, -100, -100, -100,# rep(-5, 13),
           -5, -5),
    ub = c(1000, 100, 100, 100, 100, 100, 100, 100, 100,# rep(5, 13),
           1000, 100, 100, 100, 100, 100, 100, 100, 100, 100,# rep(5, 13),
           5, 5)
  ) %>% export_all(folder = '03_Eredmenyek/masodik_proba_SA_2013_5pc')
# Ennek eredménye konvergált, de a sigma_demand a határon van! Ettől függetlenül használjuk ennek az eredményét 
# kiindulópontként az industry dummy-s modellhez, ahol ráadásul szélesítsük a sigma-k sávját!

initvec <- c(9.706202855,
             4.379259446,
             -0.048826437,
             2.742325701,
             -4.843836928,
             0.837938939,
             2.923320344,
             10.34618162,
             rep(0, 14),
             -0.595302411,
             -0.625775239,
             0.114378678,
             -15.4751159,
             0.064430476,
             0.034511131,
             -1.826905507,
             5.543018124,
             0.244696787,
             rep(0, 14),
             5,
             3.367670339
)

diseq_model_cens_nocorr_bovebb <- fitdiseq(
  demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
    x12_liqasset + x12_growtha + x12_growths + c12_own + ind_factor + 0,
  supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
    x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince + ind_factor + 0,
  data = dt_model[r <= 0.05],
  init = initvec,
  optimizer = 'SA',
  control = list('verbose' = TRUE, 'max.time' = 15*60*60),
  lb = c(-100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
         -100, -100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
         -10, -10),
  ub = c( 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
          100, 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
          10, 10)
) %>% export_all(folder = '03_Eredmenyek/masodik_proba_SA_2013_5pc_inddel')

# Ez konvergált! Az eredményében azonban az egyik szórásparaméter a határon van... 
# Konvergáltassuk újra a szórásparaméterek lazításával, de még kis mintán! 

diseq_model_cens_nocorr_bovebb <- refitdiseq(diseq_model_cens_nocorr_bovebb,
  demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
    x12_liqasset + x12_growtha + x12_growths + c12_own + ind_factor + 0,
  supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
    x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince + ind_factor + 0,
  data = dt_model[r <= 0.05],
  optimizer = 'SA',
  control = list('verbose' = TRUE, 'max.time' = 5.5*60*60),
  lb = c(-100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
         -100, -100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
         -20, -20),
  ub = c( 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
          100, 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
          20, 20)
) %>% export_all(folder = '03_Eredmenyek/masodik_proba_SA_2013_5pc_inddel_v2')

diseq_model_cens_nocorr_bovebb <- refitdiseq(diseq_model_cens_nocorr_bovebb,
                                             demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
                                               x12_liqasset + x12_growtha + x12_growths + c12_own + ind_factor + 0,
                                             supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
                                               x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince + ind_factor + 0,
                                             data = dt_model[r <= 0.1],
                                             optimizer = 'SA',
                                             control = list('verbose' = TRUE, 'max.time' = 5*60*60),
                                             lb = c(-100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
                                                    -100, -100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
                                                    -20, -20),
                                             ub = c( 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
                                                     100, 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
                                                     20, 20)
) %>% export_all(folder = '03_Eredmenyek/masodik_proba_SA_2013_10pc_inddel') %>% 
  # Nem ért el konvergenciát, ráadásul a paraméterértékek több helyen is a tartományaik határán vannak...
  # Bővített tartományokkal indítsuk újra!
  refitdiseq(lb = c(-200, -200, -200, -200, -200, -200, -200, -200, rep(-2000, 14),
                    -200, -200, -200, -200, -200, -200, -200, -200, -200, rep(-2000, 14),
                    -20, -20),
             ub = c( 200, 200, 200, 200, 200, 200, 200, 200, rep(2000, 14),
                     200, 200, 200, 200, 200, 200, 200, 200, 200, rep(2000, 14),
                     20, 20)
  ) %>% export_all(folder = '03_Eredmenyek/masodik_proba_SA_2013_10pc_inddel_2')
# Nem konvergált, és ugyan a szórásparaméterek nem tűnnek hihetetlennek, a keresleti egyenlet együtthatói elég hatalmasak, 
# úgyhogy további próbálkozásunkig haggyuk a kérdést! 

##########################################################
#######     DE módszerrel történő optimalizáció    #######
##########################################################

# Először menjen industry dummy-k nélkül a dolog. Mivel a hétvégén nem tudunk ránézni,
# nem tudjuk közben bővíteni a paraméterteret, vagy a konvergált értékeket kezdeti értékeknek használva
# áttérni a bővebb modellre. Így a szűkebb modellben maradva a mintát próbáljuk bővíteni:

diseq_model_cens_nocorr_de <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.05], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 20000, storepopfrom = 1, storepopfreq = 100),
    lb = c(-1000, -100, -100, -100, -100, -100, -100, -100, -100,# rep(-5, 13),
           -1000, -100, -100, -100, -100, -100, -100, -100, -100, -100,# rep(-5, 13),
           -10, -10),
    ub = c(1000, 100, 100, 100, 100, 100, 100, 100, 100,# rep(5, 13),
           1000, 100, 100, 100, 100, 100, 100, 100, 100, 100,# rep(5, 13),
           10, 10)
  ) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_5pc') %>% 
  refitdiseq(data = dt_model[r <= 0.1],
             control = DEoptim.control(trace = 100, itermax = 5000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_1') %>% 
  refitdiseq(data = dt_model[r <= 0.1],
             control = DEoptim.control(trace = 100, itermax = 5000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_2') %>% 
  refitdiseq(data = dt_model[r <= 0.1],
             control = DEoptim.control(trace = 100, itermax = 5000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_3') %>% 
  refitdiseq(data = dt_model[r <= 0.1],
             control = DEoptim.control(trace = 100, itermax = 5000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_4') %>% 
  refitdiseq(lb = c(-1000, -100, -100, -100, -100, -100, -100, -100, -100,# rep(-5, 13),
                    -1000, -100, -100, -100, -100, -100, -100, -100, -100, -100,# rep(-5, 13),
                    -10, -10),
             ub = c(1000, 100, 100, 100, 100, 100, 100, 100, 100,# rep(5, 13),
                    1000, 100, 100, 100, 100, 100, 100, 100, 100, 100,# rep(5, 13),
                    10, 10),
             control = DEoptim.control(trace = 100, itermax = 5000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_5') %>% 
  refitdiseq(control = DEoptim.control(trace = 100, itermax = 10000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_6') %>% 
  refitdiseq(control = DEoptim.control(trace = 100, itermax = 20000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_7')
# Nem konvergált, de alapvetően nem nagyon mászik sehova...

# Próbáljuk inkább meg ennek a kezdeti értékétől az industry dummy-s verziót!
initvec <- c(3.594407159,
             87.35529707,
             0.226204463,
             1.66818832,
             -1.79066275,
             0.247860883,
             0.919723319,
             99.99999097,
             rep(0, 14),
             -0.547598502,
             -3.664953726,
             0.151266083,
             -17.1421743,
             0.036654716,
             0.048180579,
             -0.904686941,
             5.488065493,
             0.144661622,
             rep(0, 14),
             5.286628491,
             3.305375232
)

diseq_model_cens_nocorr_bovebb_de <- fitdiseq(demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
                                               x12_liqasset + x12_growtha + x12_growths + c12_own + ind_factor + 0,
                                             supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
                                               x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince + ind_factor + 0,
                                             data = dt_model[r <= 0.1],
                                             init = initvec,
                                             optimizer = 'DE',
                                             control = DEoptim.control(trace = 100, itermax = 20000),
                                             lb = c(-100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
                                                    -100, -100, -100, -100, -100, -100, -100, -100, -100, rep(-1000, 14),
                                                    -10, -10),
                                             ub = c( 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
                                                     100, 100, 100, 100, 100, 100, 100, 100, 100, rep(1000, 14),
                                                     10, 10)
) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_inddel_1') %>% 
  # Egy napig futott, 20 ezer iteráción, és nem konvergált...
  refitdiseq(control = DEoptim.control(trace = 100, itermax = 5000)) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_inddel_2') %>% 
  refitdiseq(control = DEoptim.control(trace = 100, 
                                       itermax = 20000)) %>% 
 #                                      parallelType = 1, 
 #                                      parVar = c("y", "X", "idx_bd", "idx_bs", "idx_sd", "idx_ss")))
  export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_inddel_3')


#diseq_model_cens_nocorr_de <- readRDS('X:/PSF/_Common/Személyek/BanaiÁdám/Horváth Ákos/panel/diseq_R/03_Eredmenyek/pure_DE_2013_10pc_4/model.dem')

data <- dt_model[r <= 0.1]

diseq_model_cens_nocorr_de_10 <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.1], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 40000, storepopfrom = 1, storepopfreq = 100), #, parallelType = 1, 
                             # parVar = c("loglike.diseq.tobit")),
    lb = c(-2000, -200, -200, -200, -200, -200, -200, -200, -200,# rep(-5, 13),
           -2000, -200, -200, -200, -200, -200, -200, -200, -200, -200,# rep(-5, 13),
           -10, -10),
    ub = c(2000, 200, 200, 200, 200, 200, 200, 200, 200,# rep(5, 13),
           2000, 200, 200, 200, 200, 200, 200, 200, 200, 200,# rep(5, 13),
           10, 10)
  ) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc')

loglik_dt <- as.data.table(diseq_model_cens_nocorr_de_10$optim_trace[[1]])
loglik_dt[, iter_no := .I]
ggplot(loglik_dt) + geom_line(aes(x = iter_no, y = V1))

params_dt <- as.data.table(diseq_model_cens_nocorr_de_10$optim_trace[[2]])
params_dt[, iter_no := .I]
ggplot(params_dt) + geom_line(aes(x = iter_no, y = par1))
#alma <- refitdiseq(diseq_model_cens_nocorr_de_10)

diseq_model_cens_nocorr_de_10_2nd <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.1], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 30000, storepopfrom = 1, storepopfreq = 100), #, parallelType = 1, 
    # parVar = c("loglike.diseq.tobit")),
    lb = c(-5000, -500, -500, -500, -500, -500, -500, -500, -500,# rep(-5, 13),
           -5000, -500, -500, -500, -500, -500, -500, -500, -500, -500,# rep(-5, 13),
           0, 0),
    ub = c(5000, 500, 500, 500, 500, 500, 500, 500, 500,# rep(5, 13),
           5000, 500, 500, 500, 500, 500, 500, 500, 500, 500,# rep(5, 13),
           10, 10)
  ) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_2nd')

loglik_dt <- as.data.table(diseq_model_cens_nocorr_de_10_2nd$optim_trace[[1]])
loglik_dt[, iter_no := .I]
ggplot(loglik_dt) + geom_line(aes(x = iter_no, y = V1))

params_dt <- as.data.table(diseq_model_cens_nocorr_de_10_2nd$optim_trace[[2]])
params_dt[, iter_no := .I]
ggplot(params_dt) + geom_line(aes(x = iter_no, y = par1))


diseq_model_cens_nocorr_de_10_3rd <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.1], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 30000, storepopfrom = 1, storepopfreq = 100), #, parallelType = 1, 
    # parVar = c("loglike.diseq.tobit")),
    lb = c(-10000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,# rep(-5, 13),
           -10000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,# rep(-5, 13),
           0, 0),
    ub = c(10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,# rep(5, 13),
           10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,# rep(5, 13),
           10, 10)
  ) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_3rd')

loglik_dt <- as.data.table(diseq_model_cens_nocorr_de_10_3rd$optim_trace[[1]])
loglik_dt[, iter_no := .I]
ggplot(loglik_dt) + geom_line(aes(x = iter_no, y = V1))

params_dt <- as.data.table(diseq_model_cens_nocorr_de_10_3rd$optim_trace[[2]])
params_dt[, iter_no := .I]
ggplot(params_dt) + geom_line(aes(x = iter_no, y = par1))

# A fenti modellben a c12_own megint a paramétertér szegélyére került,
# de amíg nem jutott el oda, addig is a konstanssal együtt teljesen zűrös pályát írtak le 
# (a konstans is járt a paramétertér szélén!).
# Vegyük ki a c12_own változót, és próbáljuk újra (estére) ezzel a bővített paramétertérrel!
# Emellett a kínálati egyenletben az age és a size együtthatói gyanúsak,
# a keresleti egyenletben viszont továbbra is felfújtnak tűnnek az együtthatók!...

diseq_model_cens_nocorr_de_10_4th <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths, # + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.1], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 30000, storepopfrom = 1, storepopfreq = 100), #, parallelType = 1, 
    # parVar = c("loglike.diseq.tobit")),
    lb = c(-10000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, # -1000,# rep(-5, 13),
           -10000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,# rep(-5, 13),
           0, 0),
    ub = c(10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, # 1000,# rep(5, 13),
           10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,# rep(5, 13),
           10, 10)
  ) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_4th')

# Továbbra sem jó a dolog... Próbáljuk a következő körben a size változót is kihagyni 
# a keresleti egyenletből!!

diseq_model_cens_nocorr_de_10_5th <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths, # + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.1], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 30000, storepopfrom = 1, storepopfreq = 100), #, parallelType = 1, 
    # parVar = c("loglike.diseq.tobit")),
    lb = c(-10000, -1000, -1000, -1000, -1000, -1000, -1000, # -1000, # -1000,# rep(-5, 13),
           -10000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,# rep(-5, 13),
           0, 0),
    ub = c(10000, 1000, 1000, 1000, 1000, 1000, 1000, # 1000, # 1000,# rep(5, 13),
           10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,# rep(5, 13),
           10, 10)
  ) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_5th')

loglik_dt <- as.data.table(diseq_model_cens_nocorr_de_10_5th$optim_trace[[1]])
loglik_dt[, iter_no := .I]
ggplot(loglik_dt) + geom_line(aes(x = iter_no, y = V1))

params_dt <- as.data.table(diseq_model_cens_nocorr_de_10_5th$optim_trace[[2]])
params_dt[, iter_no := .I]
ggplot(params_dt) + geom_line(aes(x = iter_no, y = par1))

# Lehet egyet szélesíteni még a konstansok paraméterterén. Persze ettől lehet, hogy konvergen-
# ciát fogunk látni, de egyrészt nem valószínű (hiszen mielőtt elérte volna a tér szélét,
# előtte nagyon ugrált a keresleti egyenlet több együtthatója), másrészt ettől még nem lesz jó 
# a modell, hiszen az együtthatók hatalmasak és nehezen értelmezhetők!
# Próbáljuk ezért még egyszer egy szélesített paramétertérrel, c12_own nélkül, de size-zal!!!

diseq_model_cens_nocorr_de_10_6th <-
  fitdiseq(
    demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
      x12_liqasset + x12_growtha + x12_growths, # + c12_own, # + ind_factor,
    supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
      x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince, # + ind_factor,
    data = dt_model[r <= 0.1], 
    optimizer = 'DE',
    control = DEoptim.control(trace = 100, itermax = 30000, storepopfrom = 1, storepopfreq = 100), #, parallelType = 1, 
    # parVar = c("loglike.diseq.tobit")),
    lb = c(-20000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, # -1000,# rep(-5, 13),
           -20000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,# rep(-5, 13),
           0, 0),
    ub = c(20000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, # 1000,# rep(5, 13),
           20000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,# rep(5, 13),
           10, 10)
  ) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_6th')

loglik_dt <- as.data.table(diseq_model_cens_nocorr_de_10_6th$optim_trace[[1]])
loglik_dt[, iter_no := .I]
ggplot(loglik_dt) + geom_line(aes(x = iter_no, y = V1))

params_dt <- as.data.table(diseq_model_cens_nocorr_de_10_6th$optim_trace[[2]])
params_dt[, iter_no := .I]
ggplot(params_dt) + geom_line(aes(x = iter_no, y = par1))

# A helyzet továbbra is az, hogy az egyik paraméter - most a size-é - beakad, és addig nem is kezd el
# a keresleti egyenlet konvergálni...

diseq_model_cens_nocorr_bovebb_de <- fitdiseq(demand_formula = ydis13_newloans ~ x12_age + x12_size + x12_soa + x12_wcratio + 
                                                x12_liqasset + x12_growtha + x12_growths + c12_own + ind_factor + 0,
                                              supply_formula = ydis13_newloans ~ x12_age + x12_size + x12_leverage + x12_negequity +
                                                x12_roa + x12_tangasset + l12_arrears + l12_numprod + l12_banksince + ind_factor + 0,
                                              data = dt_model[r <= 0.1],
                                              optimizer = 'DE',
                                              control = DEoptim.control(trace = 100, itermax = 40000, storepopfrom = 1, storepopfreq = 100),
                                              lb = c(-200, -200, -200, -200, -200, -200, -200, -200, rep(-2000, 14),
                                                     -200, -200, -200, -200, -200, -200, -200, -200, -200, rep(-2000, 14),
                                                     -10, -10),
                                              ub = c( 200, 200, 200, 200, 200, 200, 200, 200, rep(2000, 14),
                                                      200, 200, 200, 200, 200, 200, 200, 200, 200, rep(2000, 14),
                                                      10, 10)
) %>% export_all(folder = '03_Eredmenyek/pure_DE_2013_10pc_inddel')

loglik_dt_ind <- as.data.table(diseq_model_cens_nocorr_bovebb_de$optim_trace[[1]])
loglik_dt_ind[, iter_no := .I]
ggplot(loglik_dt_ind) + geom_line(aes(x = iter_no, y = V1))

params_dt_ind <- as.data.table(diseq_model_cens_nocorr_bovebb_de$optim_trace[[2]])
params_dt_ind[, iter_no := .I]
ggplot(params_dt_ind) + geom_line(aes(x = iter_no, y = par1))


######################################
#######     Egyszerű modell    #######
######################################

diseq_model_cens_nocorr_simple <-
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
           5, 5)
  ) %>% export_all(folder = '03_Eredmenyek/simple_DE_2013_5pc')
