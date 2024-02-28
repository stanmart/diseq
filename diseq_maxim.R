library(data.table)
library(GenSA)
library(DEoptim)
library(gsubfn)
library(numDeriv)
library(mvtnorm)
library(foreach)
library(parallel)


# makes it possible to unpack lists returned by functions
list <- structure(NA, class = "result")
"[<-.result" <- function(x, ..., value) {
  args <- as.list(match.call())
  args <- args[-c(1:2, length(args))]
  length(value) <- length(args)
  for (i in seq(along = args)) {
    a <- args[[i]]
    if (!missing(a))
      eval.parent(substitute(a <- v, list(a = a, v = value[[i]])))
  }
  x
}


# Function: loglike.diseq.tobit
# Description: Calculates the log-likelihood for a disequilibrium model ala Maddala & Nelson (1974).
# Parameters:
#   - beta: Vector of coefficients.
#   - y: Vector of observed dependent variable.
#   - X: Matrix of independent variables.
#   - idx_bd: Indices of coefficients for the dependent variable mean equation.
#   - idx_bs: Indices of coefficients for the dependent variable selection equation.
#   - idx_sd: Index of the coefficient for the standard deviation of the dependent variable mean equation.
#   - idx_ss: Index of the coefficient for the standard deviation of the dependent variable selection equation.
#   - likelihood_bias: Bias term added to the likelihood to avoid taking the logarithm of zero.
# Returns: 
#   - The negative log-likelihood value.

loglike.diseq.tobit <- function (beta, y, X, idx_bd, idx_bs, idx_sd, idx_ss, likelihood_bias = 0) {
  theta1 <- X[, idx_bd, drop = FALSE] %*% beta[idx_bd]
  theta2 <- X[, idx_bs, drop = FALSE] %*% beta[idx_bs]
  sigma1 <- exp(beta[idx_sd])
  if (is.null(idx_ss)) {
    # when the sigmas are fixed to be equal
    sigma2 <- sigma1
  } else {
    sigma2 <- exp(beta[idx_ss])
  }
  f1 <- dnorm(y, mean=theta1, sd=sigma1)
  f2 <- dnorm(y, mean=theta2, sd=sigma2)
  F1 <- 1 - pnorm((y-theta1) / sigma1)
  F2 <- 1 - pnorm((y-theta2) / sigma2)
  G <- ifelse(y > 0, f1*F2 + f2*F1, 1 - F1*F2)
  -sum(log(likelihood_bias + G))
}


# Function: mvnorm_exact
# Description: Calculates the joint probabilities of two normal random variables being below zero.
# Parameters:
#   - theta1: Mean of the first normal random variable.
#   - sigma1: Standard deviation of the first normal random variable.
#   - theta2: Mean of the second normal random variable.
#   - sigma2: Standard deviation of the second normal random variable.
#   - rho: Correlation coefficient of the two normal random variables.
# Returns:
#   - The joint probabilities of the two normal random variables being below zero.

mvnorm_exact <- function(theta1, sigma1, theta2, sigma2, rho) {
  mvnorm_scalar_fun <- function(theta1, sigma1, theta2, sigma2, rho) {
    pmvnorm(lower = c(-Inf, -Inf),
            upper = c(0, 0),
            mean = c(theta1, theta2),
            sigma = matrix(c(sigma1^2, sigma1*sigma2*rho,
                             sigma1*sigma2*rho, sigma2^2),
                           nrow = 2)
    )[1]
  }

  mapply(mvnorm_scalar_fun, theta1, sigma1, theta2, sigma2, rho)
}


# Function: mvnorm_exact
# Description: Calculates the approximate joint probabilities of two normal random variables being below zero. Faster than mvnorm_exact.
# Parameters:
#   - theta1: Mean of the first normal random variable.
#   - sigma1: Standard deviation of the first normal random variable.
#   - theta2: Mean of the second normal random variable.
#   - sigma2: Standard deviation of the second normal random variable.
#   - rho: Correlation coefficient of the two normal random variables.
# Returns:
#   - The joint probabilities of the two normal random variables being below zero.

mvnorm_approx <- function(theta1, sigma1, theta2, sigma2, rho) {
  a <- -theta1 / sigma1
  b <- -theta2 / sigma2
  c <- ifelse(abs(a) >= abs(b), a, b)
  d <- ifelse(abs(a) >= abs(b), b, a)

  h <- ifelse(c <= 0, c, -c)
  k <- d
  r <- ifelse(c <= 0, rho, -rho)

  mu_approx <- -r * ifelse(h > -30, dnorm(h) / pnorm(h), -h)
  sigma_approx <- 1 + r * h * mu_approx - mu_approx ^ 2
  B <- pnorm(h) * pnorm(k, mean = mu_approx, sd = sqrt(sigma_approx))
  F_00 <- ifelse(c <= 0, B, pnorm(k) - B)

  return(F_00)
}


# Function: loglike.diseq.tobit.corr
# Description: Calculates the log-likelihood for a disequilibrium model with correlated errors.
# Parameters:
#   - beta: Vector of coefficients.
#   - y: Vector of observed dependent variable.
#   - X: Matrix of independent variables.
#   - idx_bd: Indices of coefficients for the dependent variable mean equation.
#   - idx_bs: Indices of coefficients for the dependent variable selection equation.
#   - idx_sd: Index of the coefficient for the standard deviation of the dependent variable mean equation.
#   - idx_ss: Index of the coefficient for the standard deviation of the dependent variable selection equation.
#   - likelihood_bias: Bias term added to the likelihood to avoid taking the logarithm of zero.
# Returns: 
#   - The negative log-likelihood value.

loglike.diseq.tobit.corr <- function (beta, y, X, idx_bd, idx_bs, idx_sd, idx_ss, idx_corr, likelihood_bias = 0) {
  theta1 <- X[, idx_bd, drop = FALSE] %*% beta[idx_bd]
  theta2 <- X[, idx_bs, drop = FALSE] %*% beta[idx_bs]
  sigma1 <- exp(beta[idx_sd])
  if (is.null(idx_ss)) {
    # when the sigmas are fixed to be equal
    sigma2 <- sigma1
  } else {
    sigma2 <- exp(beta[idx_ss])
  }
  rho <- tanh(beta[idx_corr])

  # y > 0 eset (F1, F2 feltételes eloszlások) - Maddala & Nelson
  f1 <- dnorm(y, mean=theta1, sd=sigma1)
  f2 <- dnorm(y, mean=theta2, sd=sigma2)
  F1 <- 1 - pnorm(y,
                  mean = theta1 + sigma1 / sigma2 * rho * (y - theta2),
                  sd = sqrt(1 - rho^2) * sigma1)
  F2 <- 1 - pnorm(y,
                  mean = theta2 + sigma2 / sigma1 * rho * (y - theta1),
                  sd = sqrt(1 - rho^2) * sigma2)

  # y < 0 esetre együttes eloszlás közelítése - Mee & Owen 1982
  F_00 <- mvnorm_approx(theta1, sigma1, theta2, sigma2, rho)
  G <- ifelse(
    y > 0,
    f1*F2 + f2*F1,
    pnorm(0, mean=theta1, sd=sigma1) + pnorm(0, mean=theta2, sd=sigma2) - F_00
  )
  -sum(log(likelihood_bias + G))
  # if (is.na(ret)) {
  #   return(Inf)
  # } else {
  #   return(ret)
  # }
  # From DEOptim docs: Note that DEoptim stops if any NA or NaN value is obtained.
  # You have to redefine your function to handle these values (for instance, set
  # NA to Inf in your objective function).

}


# Function: model.matrix.diseq
# Description: Creates the model matrix and the corresponding coefficient indices for a disequilibrium model.
# Parameters:
#   - demand_formula: Formula for the demand equation.
#   - supply_formula: Formula for the supply equation.
#   - data: Data frame containing the variables in the formulas.
#   - corr: Whether the model has correlated errors.
#   - equal_sigmas: Whether the errors on the demand and supply sides are assumed to have equal standard deviations.
# Returns:
#   - A list containing the model matrix and the coefficient indices.

model.matrix.diseq <- function (demand_formula, supply_formula, data, corr=FALSE, equal_sigmas=FALSE) {
  X_d <- model.matrix(formula(demand_formula), data = data)
  X_s <- model.matrix(formula(supply_formula), data = data)
  rows <- intersect(rownames(X_d), rownames(X_s))
  n_d <- ncol(X_d)
  n_s <- ncol(X_s)
  coef_indices <- list(
    'beta_demand' = 1 : n_d,
    'beta_supply' = (n_d+1) : (n_d+n_s),
    'sigma_demand' = n_d+n_s+1,
    'sigma_supply' = if (equal_sigmas) {NULL} else {n_d+n_s+2},
    'sigma_corr' = if (!corr) {NULL} else if (equal_sigmas) {n_d+n_s+2} else {n_d+n_s+3}
  )
  # Unavailable indices are NULL
  list(cbind(X_d[rows, ], X_s[rows, ]), coef_indices)
}


# Function: fitdiseq
# Description: Fits a disequilibrium model to the data.
# Parameters:
#   - demand_formula: Formula for the demand equation.
#   - supply_formula: Formula for the supply equation.
#   - data: Data frame containing the variables in the formulas.
#   - lb: Lower bounds for the coefficients.
#   - ub: Upper bounds for the coefficients.
#   - init: Initial values for the coefficients.
#   - initpop: Initial population for the DE algorithm.
#   - equal_sigmas: Whether the errors on the demand and supply sides are assumed to have equal standard deviations.
#   - corr: Whether the model has correlated errors.
#   - optimizer: The optimization algorithm to use. Can be 'SA', 'DE', or 'optim'.
#   - control: Control parameters for the optimization algorithm.
#   - method: The optimization method to use when optimizer is 'optim'.
#   - na.action: The NA handling function to use.
#   - random_seed: The random seed to use.
#   - elapsed_times: A list of elapsed times from previous optimization runs.
#   - prev_history: A list of optimization traces from previous optimization runs.
#   - fixed_params: A vector of fixed parameters.
#   - likelihood_bias: A bias term added to the likelihood to avoid taking the logarithm of zero.
# Returns:
#   - The fitted disequilibrium model (diseq object).

fitdiseq <- function(demand_formula,
                     supply_formula,
                     data,
                     lb = NULL,
                     ub = NULL,
                     init = NULL,
                     initpop = NULL,
                     equal_sigmas = FALSE,
                     corr = FALSE,
                     optimizer = 'SA',
                     control = (if (optimizer == 'SA') {list('verbose' = TRUE, 'max.time' = 1200)}
                                else if (optimizer == 'DE') {DEoptim.control(trace = TRUE, itermax = 1000)}
                                else if (optimizer == 'optim') {list('trace' = TRUE)}),
                     method = 'Nelder-Mead',
                     na.action = na.exclude,
                     random_seed = 1991,
                     elapsed_times = list(),
                     prev_history = list(),
                     fixed_params = NULL,
                     likelihood_bias = 0
                     ) {

  cl <- match.call()

  mf <- na.action(data[, union(all.vars(demand_formula), all.vars(supply_formula)), with=FALSE])
  attr(mf, 'demand_terms') <- terms(demand_formula)
  attr(mf, 'supply_terms') <- terms(supply_formula)
  y <- mf[[all.vars(demand_formula[[2]])]]
  list[X, coef_indices] <- model.matrix.diseq(demand_formula, supply_formula, data=mf,
                                              corr=corr, equal_sigmas=equal_sigmas)
  idx_bd <- coef_indices[['beta_demand']]
  idx_bs <- coef_indices[['beta_supply']]
  idx_sd <- coef_indices[['sigma_demand']]
  idx_ss <- coef_indices[['sigma_supply']]  # A NULL here automatically signals loglike to use sigma1 for sigma2
  idx_corr <- coef_indices[['sigma_corr']]

  num_of_params <- length(idx_bd) + length(idx_bs)
  if (equal_sigmas) {
    num_of_params <- num_of_params + 1
  } else {
     num_of_params <- num_of_params + 2
  }
  if (corr) {
    num_of_params <- num_of_params + 1
  }

  # Note: no automatic bounds will be created if not provided explicitely in the function call.
  # If the optimizer needs them, let it throw an exception.

  if (length(initpop) == 1 && initpop == FALSE) {
    # FALSE/NULL should only make a difference when calling refitdiseq
    initpop <- NULL
  }

  if (!is.null(lb) && length(lb) != num_of_params) {
    stop("Incorrect lower bound size")
  }
  if (!is.null(ub) && length(ub) != num_of_params) {
    stop("Incorrect upper bound size")
  }
  if (!is.null(init) && length(init) != num_of_params) {
    stop("Incorrect initial vector size")
  }
  if (!is.null(initpop) && ncol(initpop) != num_of_params) {
    stop('Incorrect size of vectors in the initial population.')
  }
  if (!is.null(fixed_params) && length(fixed_params) != num_of_params) {
    stop('Incorrect parameter fixing vector size.')
  }

  if (!corr) {
    loglike <- function(beta) loglike.diseq.tobit(beta, y, X, idx_bd, idx_bs, idx_sd, idx_ss, likelihood_bias)
  } else {
    loglike <- function(beta) loglike.diseq.tobit.corr(beta, y, X, idx_bd, idx_bs, idx_sd, idx_ss, idx_corr, likelihood_bias)
  }

  orig_init <- init
  orig_lb <- lb
  orig_ub <- ub
  orig_initpop <- initpop
  if (!is.null(fixed_params)) {
    if (is.null(init)) {
      stop("An initial vector is required when fixing parameters")
    }
    objective_fun <- function(beta_restr) {
      beta_full <- orig_init
      beta_full[!fixed_params] <- beta_restr
      loglike(beta_full)
    }
    # NULL indexed by anything is still NULL, so no change needed here
    init <- orig_init[!fixed_params]
    lb <- orig_lb[!fixed_params]
    ub <- orig_ub[!fixed_params]
    if (!is.null(initpop)) {
      initpop <- orig_initpop[, !fixed_params]
    }
  } else {
    objective_fun <- loglike
  }

  set.seed(random_seed)
  
  t0 <- Sys.time()
  cat(paste('Starting estimation at', t0, '\n'))

  if (optimizer == 'SA') {

    list[neg_log_likelihood_opt, beta_opt, history] <- GenSA(init, objective_fun, lb, ub, control = control)

  } else if (optimizer == 'DE') {

    if(!is.null(initpop)) {
      control$initialpop <- initpop
    } else if(!is.null(init)) {
      NP <- ifelse(!is.na(control$NP), control$NP, 10 * length(init))
      initrand <- runif((NP - 1) * length(init),
                        min = rep(lb, times = NP - 1),
                        max = rep(ub, times = NP - 1))
      initmat <- matrix(c(init, initrand),
                        ncol = length(init),
                        nrow = NP,
                        byrow = TRUE)
      control$initialpop <- initmat
      control$NP <- NP
    }

    list[opt, member] <- DEoptim(objective_fun, lb, ub, control = control)
    neg_log_likelihood_opt <- opt$bestval
    beta_opt <- opt$bestmem
    history <- list(member$bestvalit, member$bestmemit, member$pop, member$storepop)

  } else if (optimizer == 'optim') {

    if (is.null(init)) {
      cat("No initial vector provided. Using all zeros as the starting point.")
      init <- rep(0, num_of_params)
      # Therefore the starting point for sigmas is exp(0) = 1
    }
    list[beta_opt, neg_log_likelihood_opt, counts, convergence] <-
      optim(init, objective_fun, method = method, control = control)
    history <- c(counts, convergence)
    names(history) = c(names(counts), 'convergence')

  } else {

    stop('Unknown optimizer. Use one of: SA, DE, optim.')

  }
  
  t1 <- Sys.time()
  cat(paste('Estimation finished at', t1, '\n'))
  t_delta = difftime(t1, t0)

  if (!is.null(fixed_params)) {
    beta_opt_restr <- beta_opt
    beta_opt <- orig_init
    beta_opt[!fixed_params] <- beta_opt_restr
  }

  beta_names <- c(
    paste0(colnames(X)[idx_bd], ' - d'),
    paste0(colnames(X)[idx_bs], ' - s')
  )
  if (equal_sigmas) {
    beta_names <- c(beta_names, 'sigma')
  } else {
    beta_names <- c(beta_names, 'sigma_demand', 'sigma_supply')
  }
  if (corr) {
    beta_names <- c(beta_names, 'rho')
  }
  
  names(beta_opt) <- beta_names

  coefficients <- beta_opt
  coefficients[idx_sd] <- exp(coefficients[idx_sd])
  if (!is.null(idx_ss)) {
    coefficients[idx_ss] <- exp(coefficients[idx_ss])
  }
  if (corr) {
    coefficients[idx_corr] <- tanh(coefficients[idx_corr])
  }

  diseq_obj <- list(
    'coefficients' = coefficients,
    'raw_coefficients' = beta_opt,
    'coef_indices' = coef_indices,
    'call' = cl,
    'demand_terms' = terms(demand_formula),
    'supply_terms' = terms(supply_formula),
    'model' = mf,
    'log_likelihood' = -neg_log_likelihood_opt,
    'orig_rownames' = rownames(data),
    'N' = nrow(mf),
    'optim_trace' = c(prev_history, list(history)),
    'settings' = list(
      'lb' = orig_lb,
      'ub' = orig_ub,
      'init' = orig_init,
      'initpop' = orig_initpop,
      'optimizer' = optimizer,
      'control' = control,
      'method' = method,
      'corr' = corr,
      'equal_sigmas' = equal_sigmas,
      'likelihood_bias' = likelihood_bias
    ),
    'na.action' = na.action,
    'random_seed' = random_seed,
    'elapsed_times' = c(elapsed_times, list(t_delta))
  )
  attr(diseq_obj, 'class') <- 'diseq'
  return(diseq_obj)

}


# Function: refitdiseq
# Description: Refits a disequilibrium model to the data. Non-specified parameters are taken from the original model.
# Parameters:
#   - diseq_obj: The disequilibrium model to refit.
#   - demand_formula: Formula for the demand equation.
#   - supply_formula: Formula for the supply equation.
#   - data: Data frame containing the variables in the formulas.
#   - lb: Lower bounds for the coefficients.
#   - ub: Upper bounds for the coefficients.
#   - init: Initial values for the coefficients.
#   - initpop: Initial population for the DE algorithm.
#   - corr: Whether the model has correlated errors.
#   - equal_sigmas: Whether the errors on the demand and supply sides are assumed to have equal standard deviations.
#   - optimizer: The optimization algorithm to use. Can be 'SA', 'DE', or 'optim'.
#   - control: Control parameters for the optimization algorithm.
#   - method: The optimization method to use when optimizer is 'optim'.
#   - na.action: The NA handling function to use.
#   - random_seed: The random seed to use.
#   - elapsed_times: A list of elapsed times from previous optimization runs.
#   - prev_history: A list of optimization traces from previous optimization runs.
#   - fixed_params: A vector of fixed parameters.
#   - likelihood_bias: A bias term added to the likelihood to avoid taking the logarithm of zero.
#   - continue: Whether to continue from the parameters of the previous model.

refitdiseq <- function(diseq_obj,
                       demand_formula = formula(diseq_obj$demand_terms),
                       supply_formula = formula(diseq_obj$supply_terms),
                       data = diseq_obj$model,
                       lb = diseq_obj$settings$lb,
                       ub = diseq_obj$settings$ub,
                       init = NULL,
                       initpop = NULL,
                       corr = diseq_obj$settings$corr,
                       equal_sigmas = diseq_obj$settings$equal_sigmas,
                       optimizer = diseq_obj$settings$optimizer,
                       control = diseq_obj$settings$control,
                       method = diseq_obj$settings$method,
                       na.action = diseq_obj$na.action,
                       random_seed = diseq_obj$random_seed,
                       elapsed_times = NULL,
                       prev_history = NULL,
                       fixed_params = NULL,
                       likelihood_bias = diseq_obj$settings$likelihood_bias,
                       continue = TRUE
                       ) {

  if (is.null(init)) {
    if (continue) {
      init <- diseq_obj$raw_coefficients
    } else {
      init <- diseq_obj$settings$init
    }
  }

  if (optimizer == "DE" && is.null(initpop)) {
    if (continue) {
      n <- length(diseq_obj$optim_trace)
      if (length(initpop) == 1 && initpop == FALSE) {
        initpop <- NULL
      } else if (!is.null(diseq_obj$optim_trace[[n]][[3]])) {
        initpop <- initpop <- diseq_obj$optim_trace[[n]][[3]]
      }
    } else {
      initpop <- diseq_obj$settings$initpop
    }
  }

  if (is.null(elapsed_times)) {
    if (continue) {
      # Append the elapsed times with a new value
      elapsed_times <- diseq_obj$elapsed_times
    } else {
      # Replace the last elapsed time with a new value
      elapsed_times <- diseq_obj$elapsed_times[-length(diseq_obj$elapsed_times)]
    }
  }

  if (is.null(prev_history)) {
    if (continue) {
      # Append the history with a new value
      prev_history <- diseq_obj$optim_trace
    } else {
      # Replace the last step of the history with a new value
      prev_history <- diseq_obj$optim_trace[-length(diseq_obj$elapsed_times)]
    }
  }

  # The next three ifs are necessary for backward compatibility reasons
  if(is.null(corr)) {
    corr <- FALSE
  }
  if (is.null(equal_sigmas)) {
    equal_sigmas <- FALSE
  }
  if(is.null(likelihood_bias)) {
    likelihood_bias <- 0
  }

  fitdiseq(
    demand_formula = demand_formula,
    supply_formula = supply_formula,
    data = data,
    lb = lb,
    ub = ub,
    init = init,
    initpop = initpop,
    corr = corr,
    equal_sigmas = equal_sigmas,
    optimizer = optimizer,
    control = control,
    method = method,
    na.action = na.action,
    random_seed = random_seed,
    elapsed_times = elapsed_times,
    prev_history = prev_history,
    fixed_params = fixed_params,
    likelihood_bias = likelihood_bias
  )

}


# Function: ll_contributions
# Description: Calculates the contribution of each observation to the log-likelihood.
# Parameters:
#   - demand_formula: Formula for the demand equation.
#   - supply_formula: Formula for the supply equation.
#   - data: Data frame containing the variables in the formulas.
#   - beta: Vector of coefficients.
#   - equal_sigmas: Whether the errors on the demand and supply sides are assumed to have equal standard deviations.
#   - na.action: The NA handling function to use.
#   - corr: Whether the model has correlated errors.
#   - likelihood_bias: A bias term added to the likelihood to avoid taking the logarithm of zero.
# Returns:
#   - A data table containing the contributions to the log-likelihood.

ll_contributions <- function(demand_formula = NULL,
                             supply_formula = NULL,
                             data = NULL,
                             beta = NULL,
                             equal_sigmas = FALSE,                             
                             na.action = na.exclude,
                             corr = FALSE,
                             likelihood_bias = 0
                             ) {

  mf <- na.action(data[, union(all.vars(demand_formula), all.vars(supply_formula)), with=FALSE])
  attr(mf, 'demand_terms') <- terms(demand_formula)
  attr(mf, 'supply_terms') <- terms(supply_formula)
  y <- mf[[all.vars(demand_formula[[2]])]]
  list[X, coef_indices] <- model.matrix.diseq(demand_formula, supply_formula, data=mf,
                                              corr=corr, equal_sigmas=equal_sigmas)
  idx_bd <- coef_indices[['beta_demand']]
  idx_bs <- coef_indices[['beta_supply']]
  idx_sd <- coef_indices[['sigma_demand']]
  idx_ss <- coef_indices[['sigma_supply']]
  idx_corr <- coef_indices[['sigma_corr']]

  y_name <- as.character(demand_formula[[2]])
  
  contributions <- as.data.table(X %*% diag(beta[c(idx_bd, idx_bs)]))
  setnames(contributions, c(paste0('beta * ', colnames(X)[idx_bd], ' - d'),
                            paste0('beta * ', colnames(X)[idx_bs], ' - s')))
  contributions[, (y_name) := y]

  theta1 <- X[, idx_bd] %*% beta[idx_bd]
  theta2 <- X[, idx_bs] %*% beta[idx_bs]
  sigma1 <- exp(beta[idx_sd])
  if (is.null(idx_ss)) {
    sigma2 <- sigma1
  } else {
    sigma2 <- exp(beta[idx_ss])
  }
  
  if (corr) {
    rho <- tanh(beta[idx_corr])
  }

  if (!corr) {
    f1 <- dnorm(y, mean=theta1, sd=sigma1)
    f2 <- dnorm(y, mean=theta2, sd=sigma2)
    F1 <- 1 - pnorm((y-theta1) / sigma1)
    F2 <- 1 - pnorm((y-theta2) / sigma2)
    prob_zero <- 1 - F1*F2
    G <- ifelse(y > 0, f1*F2 + f2*F1, prob_zero)
  } else {
    f1 <- dnorm(y, mean=theta1, sd=sigma1)
    f2 <- dnorm(y, mean=theta2, sd=sigma2)
    F1 <- 1 - pnorm(y,
                    mean = theta1 + sigma1 / sigma2 * rho * (y - theta2),
                    sd = sqrt(1 - rho^2) * sigma1)
    F2 <- 1 - pnorm(y,
                    mean = theta2 + sigma2 / sigma1 * rho * (y - theta1),
                    sd = sqrt(1 - rho^2) * sigma2)

    # y < 0 esetre együttes eloszlás közelítése - Mee & Owen 1982
    F_00 <- mvnorm_approx(theta1, sigma1, theta2, sigma2, rho)
    prob_zero <- pnorm(0, mean=theta1, sd=sigma1) + pnorm(0, mean=theta2, sd=sigma2) - F_00
    G <- ifelse(
      y > 0,
      f1*F2 + f2*F1,
      prob_zero
    )
  }

  contributions[, pred_demand := theta1]
  contributions[, pred_supply := theta2]
  contributions[, pred_outcome := pmin(pred_demand, pred_supply)]

  contributions[, `demand_density (f1)` := f1]
  contributions[, `supply_density (f2)` := f2]
  contributions[, `prob_demand_above_observed (F1)` := F1]
  contributions[, `prob_supply_above_observed (F2)` := F2]
  contributions[, prob_zero := prob_zero]

  contributions[, likelihood_contribution := G]
  contributions[, loglike_contribution := log(G)]
  contributions[, biased_loglike_contribution := log(G + likelihood_bias)]

  return(contributions)

}


# Function: model_ll_contributions
# Description: Convenience function to calculate the log-likelihood contributions for a fitted disequilibrium model.
# Parameters:
#   - diseq_obj: The fitted disequilibrium model.
# Returns:
#   - A data table containing the contributions to the log-likelihood.

model_ll_contributions <- function(diseq_obj) {

  corr <- !is.null(diseq_obj$settings$corr) && diseq_obj$settings$corr == TRUE
  equal_sigmas <- !is.null(diseq_obj$settings$equal_sigmas) && diseq_obj$settings$equal_sigmas == TRUE

  ll_contributions(demand_formula = formula(diseq_obj$demand_terms),
                   supply_formula = formula(diseq_obj$supply_terms),
                   data = diseq_obj$model,
                   beta = diseq_obj$raw_coefficients,
                   corr = corr,
                   equal_sigmas = equal_sigmas,
                   na.action = diseq_obj$na.action,
                   likelihood_bias = diseq_obj$settings$likelihood_bias
                   )

}


# Function: print.diseq
# Description: Prints the fitted disequilibrium model, including the demand and supply equations and the coefficients.

print.diseq <- function(diseq_obj) {

  cat('\nCall:\n')
  print(diseq_obj$call)
  cat('\n')

  cat('Demand equation:\n')
  print(diseq_obj$coefficients[ diseq_obj$coef_indices[['beta_demand']] ])
  print(diseq_obj$coefficients[ diseq_obj$coef_indices[['sigma_demand']] ])
  cat('\n')

  cat('Supply equation:\n')
  print(diseq_obj$coefficients[ diseq_obj$coef_indices[['beta_supply']] ])
  if (diseq_obj$settings$equal_sigmas) {
    print(diseq_obj$coefficients[ diseq_obj$coef_indices[['sigma_demand']] ])
  } else {
    print(diseq_obj$coefficients[ diseq_obj$coef_indices[['sigma_supply']] ])
  }
  cat('\n')

  if (!is.null(diseq_obj$settings$corr) && diseq_obj$settings$corr) {
    cat('Correlation of error terms:\n')
    print(diseq_obj$coefficients[ diseq_obj$coef_indices[['sigma_corr']] ])
    cat('\n')
  }

}


# Function: predict.diseq
# Description: Predicts from a fitted disequilibrium model.
# Parameters:
#   - diseq_obj: The fitted disequilibrium model.
#   - newdata: Data frame containing the variables in the formulas. If NULL, the original data is used.
#   - type: The type of prediction to make. Can be 'prob_supply_constrained', 'demand', 'supply', 'response', or 'expected_deficit'.
#   - conditional: Whether to make conditional predictions.
#   - na.action: The NA handling function to use.
#   - prob_without_positive_demand: Whether negative demand values should be included in the probability of supply constraint. Deprecated.
#   - exact: Whether to use the exact or approximate method for the joint probabilities.
# Returns:
#   - The predicted values.

predict.diseq <- function(diseq_obj,
                          newdata = NULL,
                          type='prob_supply_constrained',
                          conditional=FALSE,
                          na.action = diseq_obj$na.action,
                          prob_without_positive_demand = FALSE,
                          exact = TRUE) {

  corr <- !is.null(diseq_obj$settings$corr) && diseq_obj$settings$corr == TRUE
  equal_sigmas <- !is.null(diseq_obj$settings$equal_sigmas) && diseq_obj$settings$equal_sigmas == TRUE

  if (is.null(newdata)) {

    list[X, coef_indices] = model.matrix.diseq(
      diseq_obj$demand_terms,
      diseq_obj$supply_terms,
      data = diseq_obj$model,
      corr = corr,
      equal_sigmas = equal_sigmas
    )
    if (conditional) {
      y <- diseq_obj$model[[all.vars(formula(diseq_obj$demand_terms)[[2]])]] 
    }
    orig_rownames <- diseq_obj$orig_rownames

  } else {

    if (!conditional) {
      vars_to_select <- union(all.vars(update(diseq_obj$demand_terms, NULL ~ .)),
                              all.vars(update(diseq_obj$supply_terms, NULL ~ .)))
    } else {
      vars_to_select <- union(all.vars(diseq_obj$demand_terms),
                              all.vars(diseq_obj$supply_terms))
    }

    newdata_na_cleaned <- na.action(newdata[, vars_to_select, with=FALSE])
    list[X, coef_indices] = model.matrix.diseq(
      update(diseq_obj$demand_terms, NULL ~ .),
      update(diseq_obj$supply_terms, NULL ~ .),
      data = newdata_na_cleaned,
      corr = corr,
      equal_sigmas = equal_sigmas
    )
    if (conditional) {
      y <- newdata_na_cleaned[[all.vars(formula(diseq_obj$demand_terms)[[2]])]] 
    }
    orig_rownames <- rownames(newdata)

  }
  new_rownames <- rownames(X)

  idx_bd <- coef_indices[['beta_demand']]
  idx_bs <- coef_indices[['beta_supply']]
  idx_sd <- coef_indices[['sigma_demand']]
  idx_ss <- coef_indices[['sigma_supply']]
  idx_corr <- coef_indices[['sigma_corr']]

  beta_opt <- diseq_obj$raw_coefficients

  theta_demand <- X[, idx_bd] %*% beta_opt[idx_bd]
  theta_supply <- X[, idx_bs] %*% beta_opt[idx_bs]

  sigma1 <- exp(beta_opt[idx_sd])
  if (is.null(idx_ss)) {
    sigma2 <- sigma1
  } else {
    sigma2 <- exp(beta_opt[idx_ss])
  }
  if (is.null(idx_corr)) {
    rho <- 0
  } else {
    rho <- tanh(beta_opt[idx_corr])
  }

  if (type == 'prob_supply_constrained' && conditional == FALSE) {
    out <- predict_prob_supply_constrained_uncond(
      theta_demand, theta_supply, sigma1, sigma2, rho,
      prob_without_positive_demand, exact
    )
  } else if (type == 'demand' && conditional == FALSE) {
    out <- theta_demand
  } else if (type == 'supply' && conditional == FALSE) {
    out <- theta_supply
  } else if (type == 'response' && conditional == FALSE) {
    out <- pmax(0, pmin(theta_supply, theta_demand))
  } else if (type == 'expected_deficit' && conditional == FALSE) {
    out <- predict_expected_deficit_uncond(
      theta_demand, theta_supply, sigma1, sigma2, rho, exact
    )
  } else if (type == 'prob_supply_constrained' & conditional == TRUE) {
    if (prob_without_positive_demand) {
      stop("prob_without_positive_demand is a legacy argument and is not implemented for conditional probabilities")
    }
    out <- predict_prob_supply_constrained_cond(
        y, theta_demand, theta_supply, sigma1, sigma2, rho, exact
    )
  } else if (type == 'expected_deficit' && conditional == TRUE) {
    out <- predict_expected_deficit_cond(
      y, theta_demand, theta_supply, sigma1, sigma2, rho, exact
    )
  } else if (type == 'response' && conditional == TRUE) {
    warning("Chosing type=response and conditional=TRUE returns the observed outcome")
    out <- y
  } else {
    stop('Supplied prediction type and conditional/unconditional combination is not implemented.')
  }

  if (identical(na.action, na.exclude)) {
    out_final <- rep(NA, times = length(orig_rownames))
    names(out_final) <- orig_rownames
    out_final[new_rownames] <- out
  }
  else {
    out_final <- out
  }
  return(out_final)

}


# Function: predict_prob_supply_constrained_uncond
# Description: Predicts the unconditional probability of being supply constrained.

predict_prob_supply_constrained_uncond <- function(theta_demand, theta_supply, sigma1, sigma2, rho,
                                                   prob_without_positive_demand = FALSE, exact = TRUE) {

    if (prob_without_positive_demand == TRUE) {
      # Legacy code, should not be used
      sigma <- sqrt(sigma1^2 + sigma2^2 - 2 * rho * sigma1 * sigma2)
      prob_supply_constrained <- pnorm((theta_demand-theta_supply) / sigma)

    } else {

      mvnorm_fun <- if (exact == TRUE) {mvnorm_exact} else {mvnorm_approx}

      # A P(D > S metszet D > 0) = P(S - D < 0 metszet -D < 0) együttes eloszlása
      # levezetés papíron Peti mappájában (térkép!)
      mu_1 <- theta_supply - theta_demand
      mu_2 <- -theta_demand
      sigma_1_bar <- sqrt(sigma1^2 + sigma2^2 - 2 * rho * sigma1 * sigma2)
      sigma_2_bar <- sigma2
      corr_coef <- sigma1 / sigma2 - rho

      prob_supply_constrained <- mvnorm_fun(mu_1, sigma_1_bar, mu_2, sigma_2_bar, corr_coef)

    }

    return(prob_supply_constrained)

}


# Function: predict_prob_supply_constrained_cond
# Description: Predicts the conditional probability of being supply constrained.

predict_prob_supply_constrained_cond <- function(y, theta_demand, theta_supply, sigma_d, sigma_s, rho, exact = TRUE) {

  mvnorm_fun <- if (exact == TRUE) {mvnorm_exact} else {mvnorm_approx}

  # y > 0 part
  mu_D_cond_S_y <- theta_demand + rho * sigma_d / sigma_s * (y - theta_supply)
  sigma_D_cond_S_y <- sigma_d * sqrt(1 - rho^2)

  mu_S_cond_D_y <- theta_supply + rho * sigma_s / sigma_d * (y - theta_demand)
  sigma_S_cond_D_y <- sigma_s * sqrt(1 - rho^2)

  P_y_smaller_D_cond_S_eq_y <- 1 - pnorm(y, mean = mu_D_cond_S_y, sd = sigma_D_cond_S_y)
  P_y_smaller_S_cond_D_eq_y <- 1 - pnorm(y, mean = mu_S_cond_D_y, sd = sigma_S_cond_D_y)
  f_D <- dnorm(y, mean = theta_demand, sd = sigma_d)
  f_S <- dnorm(y, mean = theta_supply, sd = sigma_s)

  prob_cond_pos_y <- (P_y_smaller_D_cond_S_eq_y * f_S) / (P_y_smaller_D_cond_S_eq_y * f_S + P_y_smaller_S_cond_D_eq_y * f_D)

  # y == 0 part
  P_D_smaller_0 <- pnorm(0, theta_demand, sigma_d)
  P_S_smaller_0 <- pnorm(0, theta_supply, sigma_s)
  P_both_smaller_0 <- mvnorm_fun(theta_demand, sigma_d, theta_supply, sigma_s, rho)

  prob_cond_zero_y <- (P_S_smaller_0 - P_both_smaller_0) / (P_S_smaller_0 + P_D_smaller_0 - P_both_smaller_0)

  prob_supply_constrained <- ifelse(y > 0, prob_cond_pos_y, prob_cond_zero_y)

  return(prob_supply_constrained)

}


# Function: predict_prob_supply_constrained_uncond
# Description: Predicts the unconditional expected deficit (demand - supply).

predict_expected_deficit_uncond <- function(theta_demand, theta_supply, sigma_d, sigma_s, rho, exact = TRUE) {

  # Preallocate result vectors
  N <- length(theta_demand)
  expected_deficit_part_below_0 <- rep(NA, times = N)
  expected_deficit_part_above_0 <- rep(NA, times = N)

  # Note: D_star = D | S = eta ~ N(m, s)

  s <- sigma_d * sqrt(1 - rho^2)

  for (i in 1 : N) {

    fun_to_integrate_below_0 <- function(eta) {
      m <- theta_demand[i] + rho * sigma_d / sigma_s * (eta - theta_supply[i])
      E_D_star_cond_D_star_greater_0 <- m + s * dnorm(m / s) / pnorm(m / s)
      P_D_star_greater_0 <- 1 - pnorm(0, mean = m, sd = s)
      f_S_eta <- dnorm(eta, mean = theta_supply[i], sd = sigma_s)
      # The ifelse below is needed to get rid of 0 * Inf issues
      # The limit is of the product is 0 at -Inf so it should be fine
      return(ifelse(
        P_D_star_greater_0 * f_S_eta == 0,
        0,
        E_D_star_cond_D_star_greater_0 * P_D_star_greater_0 * f_S_eta
      ))
    }

    fun_to_integrate_above_0 <- function(eta) {
      m <- theta_demand[i] + rho * sigma_d / sigma_s * (eta - theta_supply[i])
      E_D_star_cond_D_star_greater_eta <- m + s * dnorm((eta - m) / s) / (1 - pnorm((eta - m) / s))
      P_D_star_greater_eta <- 1 - pnorm(eta, mean = m, sd = s)
      f_S_eta <- dnorm(eta, mean = theta_supply[i], sd = sigma_s)
      # The ifelse below is needed to get rid of 0 * Inf issues
      # The limit is of the product is 0 at Inf so it should be fine
      return(ifelse(
        P_D_star_greater_eta * f_S_eta == 0,
        0,
        (E_D_star_cond_D_star_greater_eta - eta) * P_D_star_greater_eta * f_S_eta
      ))
    }

    expected_deficit_part_below_0[i] <- integrate(fun_to_integrate_below_0, -Inf, 0)$value
    expected_deficit_part_above_0[i] <- integrate(fun_to_integrate_above_0, 0, Inf)$value
  
  }

  return(expected_deficit_part_below_0 + expected_deficit_part_above_0)

}


# Function: predict_prob_supply_constrained_uncond
# Description: Predicts the conditional expected deficit (demand - supply).

predict_expected_deficit_cond <- function(y, theta_demand, theta_supply, sigma_d, sigma_s, rho, exact = TRUE) {

  prob_supply_constrained <- predict_prob_supply_constrained_cond(
    y, theta_demand, theta_supply, sigma_d, sigma_s, rho, exact = exact
  )

  # y > 0 part
  mu_D_cond_S_y <- theta_demand + rho * sigma_d / sigma_s * (y - theta_supply)
  sigma_D_cond_S_y <- sigma_d * sqrt(1 - rho^2)

  y_bar <- (y - mu_D_cond_S_y) / sigma_D_cond_S_y

  expected_deficit_cond_constr <- mu_D_cond_S_y + sigma_D_cond_S_y * dnorm(y_bar) / (1 - pnorm(y_bar))

  # y == 0 part
  s <- sigma_d * sqrt(1 - rho^2)
  for (i in seq_along(y)) {
    if (y[i] == 0) {
      fun_to_integrate <- function(z) {
        m <- theta_demand[i] + rho * sigma_d / sigma_s * (z - theta_supply[i])
        E_D_star_cond_D_star_greater_0 <- m + s * dnorm(m / s) / pnorm(m / s)
        f_S_z <- dnorm(z, mean = theta_supply[i], sd = sigma_s)
        # The ifelse below is needed to get rid of 0 * Inf issues
        # The limit is of the product is 0 at -Inf so it should be fine
        return(ifelse(f_S_z == 0, 0, E_D_star_cond_D_star_greater_0 * f_S_z))
      }
      numerator <- integrate(fun_to_integrate, -Inf, 0)$value
      denominator <- pnorm(0, theta_supply[i], sigma_s)
      expected_deficit_cond_constr[i] <- numerator / denominator
    }

  }

  return(expected_deficit_cond_constr * prob_supply_constrained)

}


# Function: score.diseq.tobit
# Description: Calculates the score function for a fitted disequilibrium model.
# Parameters:
#   - beta: The vector of coefficients.
#   - y: The observed outcome.
#   - X: The model matrix.
#   - idx_bd: Indices of the demand coefficients.
#   - idx_bs: Indices of the supply coefficients.
#   - idx_sd: Index of the demand sigma.
#   - idx_ss: Index of the supply sigma.
# Returns:
#   - The scores.

score.diseq.tobit <- function(beta, y, X, idx_bd, idx_bs, idx_sd, idx_ss) {
  theta1 <-X[, idx_bd, drop = FALSE] %*% beta[idx_bd]
  theta2 <-X[, idx_bs, drop = FALSE] %*% beta[idx_bs]
  sigma1 <- exp(beta[idx_sd])
  if (is.null(idx_ss)) {
    sigma2 <- sigma1 
  } else {
    sigma2 <- exp(beta[idx_ss])
  }

  f1 <- dnorm(y, mean=theta1, sd=sigma1)
  f2 <- dnorm(y, mean=theta2, sd=sigma2)
  F1 <- 1 - pnorm((y-theta1) / sigma1)
  F2 <- 1 - pnorm((y-theta2) / sigma2)
  G <- ifelse(y > 0, f1*F2 + f2*F1, 1 - F1*F2)
  h1 <- (y - theta1) / sigma1
  h2 <- (y - theta2) / sigma2
  
  g_bd <- ifelse(y > 0,
                 (f1*h1*F2/sigma1 + f1*f2) / G,
                 -(f1*F2) / G
                 )
  g_bs <- ifelse(y > 0,
                 (f2*h2*F1/sigma2 + f1*f2) / G,
                 -(f2*F1) / G
                 )
  g_sd <- ifelse(y > 0,
                 (f1*F2*(h1^2-1) / (2*sigma1^2) + f1*f2*h1 / (2*sigma1)) / G,
                 -(f1*h1*F2 / (2*sigma1)) / G
                 )
  g_ss <- ifelse(y > 0,
                 (f2*F1*(h2^2-1) / (2*sigma2^2) + f1*f2*h2 / (2*sigma2)) / G,
                 -(f2*h2*F1 / (2*sigma2)) / G
                 )
  
  if (is.null(idx_ss)) {
    return(cbind(
      X[, idx_bd, drop = FALSE] * g_bd,  # Chain rule for betas
      X[, idx_bs, drop = FALSE] * g_bs,  # and using the fact that vectors are replicated columnwise
      (g_sd + g_ss) * 2 * sigma1 ^2  # Chain rule for sigma
    ))
  } else {
    return(cbind(
      X[, idx_bd, drop = FALSE] * g_bd,  # Chain rule for the betas
      X[, idx_bs, drop = FALSE] * g_bs,  # and using the fact that vectors are replicated columnwise
      g_sd * 2 * sigma1 ^2,  # Chain rule for the sigmas
      g_ss * 2 * sigma2 ^2
    ))
  }
}


# Function: summary.diseq
# Description: Calculate standard errors and p values for a fitted disequilibrium model.
# Parameters:
#   - diseq_obj: The fitted disequilibrium model.
#   - se_type: The type of standard errors to calculate. Can be 'IM', 'score', or 'robust'.
# Returns:
#   - A summary object.

summary.diseq <- function(diseq_obj, se_type = "IM") {

  corr <- !is.null(diseq_obj$settings$corr) && diseq_obj$settings$corr == TRUE
  equal_sigmas <- !is.null(diseq_obj$settings$equal_sigmas) && diseq_obj$settings$equal_sigmas == TRUE

  y <- diseq_obj$model[[all.vars(formula(diseq_obj$demand_terms)[[2]])]]
  list[X, coef_indices] = model.matrix.diseq(
    diseq_obj$demand_terms,
    diseq_obj$supply_terms,
    data = diseq_obj$model,
    corr = corr,
    equal_sigmas = equal_sigmas
  )
  idx_bd <- coef_indices[['beta_demand']]
  idx_bs <- coef_indices[['beta_supply']]
  idx_sd <- coef_indices[['sigma_demand']]
  idx_ss <- coef_indices[['sigma_supply']]
  idx_corr <- coef_indices[['sigma_corr']]

  beta_opt <- diseq_obj$raw_coefficients

  if (!corr) {
    loglike <- function(beta) loglike.diseq.tobit(beta, y, X, idx_bd, idx_bs, idx_sd, idx_ss)
  } else {
    loglike <- function(beta) loglike.diseq.tobit.corr(beta, y, X, idx_bd, idx_bs, idx_sd, idx_ss, idx_corr)
  }

  if (se_type == "IM") {
    I_opt <- hessian(loglike, beta_opt)
    raw_std_err <- sqrt(diag(solve(I_opt)))
  } else if (se_type == "score") {
    if (corr) stop ("To be implemented in the correlated case")
    scores <- score.diseq.tobit(beta_opt, y, X, idx_bd, idx_bs, idx_sd, idx_ss)
    raw_std_err <- sqrt(diag(solve(t(scores) %*% scores)))
  } else if (se_type == "robust") {
    if (corr) stop ("To be implemented in the correlated case")
    I_opt <- hessian(loglike, beta_opt)
    scores <- score.diseq.tobit(beta_opt, y, X, idx_bd, idx_bs, idx_sd, idx_ss)
    raw_std_err <- sqrt(diag(solve(I_opt) %*% (t(scores) %*% scores) %*% solve(I_opt)))
  } else {
    stop("Unknown standard error type. Must be one of: 'IM', 'score', 'robust'")
  }

  raw_t_value <- beta_opt / raw_std_err
  raw_p_value <- pnorm(-abs(raw_t_value)) * 2

  coefficients <- diseq_obj$coefficients
  std_err <- raw_std_err
  std_err[idx_sd] <- std_err[idx_sd] * exp(beta_opt[idx_sd])  # Delta-method
  if (!is.null(idx_ss)) {
    std_err[idx_ss] <- std_err[idx_ss] * exp(beta_opt[idx_ss])  # Delta-method
  }
  if (!is.null(idx_corr)) {
    std_err[idx_corr] <- std_err[idx_corr] * (1 - tanh(beta_opt[idx_ss])^2)  # Delta-method
  }
  t_value <- coefficients / std_err
  p_value <- pnorm(-abs(t_value)) * 2

  raw_coefficients <- as.matrix(data.frame(
    'Estimate' = beta_opt,
    'Std. Error' = raw_std_err,
    't value' = raw_t_value,
    'Pr(>[t])' = raw_p_value,
    row.names = names(beta_opt)
    ))

  transformed_coefficients <- as.matrix(data.frame(
    'Estimate' = coefficients,
    'Std. Error' = std_err,
    't value' = t_value,
    'Pr(>[t])' = p_value,
    row.names = names(coefficients)
  ))

  diseq_summary_obj <- list(
    'raw_coefficients' = raw_coefficients,
    'coefficients' = transformed_coefficients,
    'call' = diseq_obj$call,
    'demand_terms' = diseq_obj$demand_terms,
    'supply_terms' = diseq_obj$supply_terms,
    'corr' = corr,
    'equal_sigmas' = equal_sigmas,
    'coef_indices' = diseq_obj$coef_indices,
    'log_likelihood' = diseq_obj$log_likelihood,
    'N' = diseq_obj$N,
    'orig_rownames' = diseq_obj$orig_rownames,
    'optim_trace' = diseq_obj$optim_trace,
    'na.action' = diseq_obj$na.action
  )
  attr(diseq_summary_obj, 'class') <- 'summary.diseq'
  return(diseq_summary_obj)

}


# Function: print.summary.diseq
# Description: Pretty-prints a summary of a fitted disequilibrium model.
# Parameters:
#   - diseq_summary_obj: The summary object.

print.summary.diseq <- function(diseq_summary_obj) {

  corr <- !is.null(diseq_summary_obj$corr) && diseq_summary_obj$corr == TRUE

  idx_bd <- diseq_summary_obj$coef_indices[['beta_demand']]
  idx_bs <- diseq_summary_obj$coef_indices[['beta_supply']]
  idx_sd <- diseq_summary_obj$coef_indices[['sigma_demand']]
  idx_ss <- diseq_summary_obj$coef_indices[['sigma_supply']]
  idx_corr <- diseq_summary_obj$coef_indices[['sigma_corr']]


  # Not very nice but works for the purposes of this function
  if (is.null(idx_ss)) {
    idx_ss <- idx_sd
  }

  est_table_demand <- as.data.frame(diseq_summary_obj$coefficients[c(idx_bd, idx_sd), ])
  est_table_supply <- as.data.frame(diseq_summary_obj$coefficients[c(idx_bs, idx_ss), ])
  if (corr) {
    est_table_corr <- as.data.frame(diseq_summary_obj$coefficients[idx_corr, , drop = FALSE])
  }

  cat('\nCall:\n')
  print(diseq_summary_obj$call)
  cat('\n')

  cat('Demand equation:\n')
  print(est_table_demand)
  cat('\n')

  cat('Supply equation:\n')
  print(est_table_supply)
  cat('\n')

  if (corr) {
    cat('Correlation of error terms:\n')
    print(est_table_corr)
    cat('\n')
  }

  cat('Number of observations: ')
  cat(diseq_summary_obj$N)
  cat('\n')
  cat('Log-likelihood: ')
  cat(diseq_summary_obj$log_likelihood)
  cat('\n\n')
}


# Function: export_to_csv
# Description: Exports a summary of a fitted disequilibrium model to a CSV file.
# Parameters:
#   - diseq_summary_obj: The summary object.
#   - file: The file path to export to.

export_to_csv <- function (diseq_summary_obj, file) {

  corr <- !is.null(diseq_summary_obj$corr) && diseq_summary_obj$corr == TRUE

  idx_bd <- diseq_summary_obj$coef_indices[['beta_demand']]
  idx_bs <- diseq_summary_obj$coef_indices[['beta_supply']]
  idx_sd <- diseq_summary_obj$coef_indices[['sigma_demand']]
  idx_ss <- diseq_summary_obj$coef_indices[['sigma_supply']]
  idx_corr <- diseq_summary_obj$coef_indices[['sigma_corr']]

  # Not very nice but works for the purposes of this function
  if (is.null(idx_ss)) {
    idx_ss <- idx_sd
  }

  est_table_demand <- as.data.frame(diseq_summary_obj$coefficients[c(idx_bd, idx_sd), ])
  est_table_supply <- as.data.frame(diseq_summary_obj$coefficients[c(idx_bs, idx_ss), ])
  if (corr) {
    est_table_corr <- as.data.frame(diseq_summary_obj$coefficients[idx_corr, , drop = FALSE])
  }

  write(
    paste('function;', toString(diseq_summary_obj$call[[1]]), ';;;', sep = ''),
    file = file, append = FALSE
  )
  for (arg in names(diseq_summary_obj$call)[-1]) {
    callstr <- tryCatch(
      paste(gsub('    ', '', format(diseq_summary_obj$call[[arg]])), collapse = ''),
      error = function (e) toString(diseq_summary_obj$call[[arg]])
    )
    write(
      paste(arg, ';', callstr, ';;;', sep = ''),
      file = file, append = TRUE)
  }

  write(';;;;', file = file, append = TRUE)
  write('Demand equation;;;;', file = file, append = TRUE)
  write(';Estimate;Std. error;t-value;p-value', file = file, append = TRUE)
  suppressWarnings(
    write.table(est_table_demand, row.names = TRUE, col.names = FALSE, file = file, append = TRUE, sep=';', na = '')
    )
  write(';;;;', file = file, append = TRUE)

  write('Supply equation;;;;', file = file, append = TRUE)
  write(';Estimate;Std. error;t-value;p-value', file = file, append = TRUE)
  suppressWarnings(
    write.table(est_table_supply, row.names = TRUE, col.names = FALSE, file = file, append = TRUE, sep=';', na = '')
    )
  write(';;;;', file = file, append = TRUE)

  if (corr) {
    write('Correlation of error terms;;;;', file = file, append = TRUE)
    write(';Estimate;Std. error;t-value;p-value', file = file, append = TRUE)
    suppressWarnings(
      write.table(est_table_corr, row.names = TRUE, col.names = FALSE, file = file, append = TRUE, sep=';', na = '')
    )
    write(';;;;', file = file, append = TRUE)
  }

  write(paste('Number of observations;', diseq_summary_obj$N, ';;;'), file = file, append = TRUE)
  write(paste('Log-likelihood;', diseq_summary_obj$log_likelihood, ';;;'), file = file, append = TRUE)
}


# Function: export_all
# Description: Exports a summary and the model object itself for a fitted disequilibrium model to a folder.
# Parameters:
#   - diseq_obj: The fitted disequilibrium model.
#   - summary_obj: The summary object. If NULL, it will be calculated.
#   - folder: The folder to export to.

export_all <- function(diseq_obj,
                       summary_obj = summary(diseq_obj),
                       folder = stop('A folder is required.')) {
  dir.create(folder)
  if (is.null(summary_obj)) {
    try({
      summary_obj <- summary(diseq_obj)
      export_to_csv(summary_obj, file.path(folder, 'summary.csv'))
    })
  } else {
    export_to_csv(summary_obj, file.path(folder, 'summary.csv'))
  }
  saveRDS(diseq_obj, file.path(folder, 'model.dem'))
  return(invisible(diseq_obj))
}

