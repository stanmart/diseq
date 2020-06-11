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
  rho <- beta[idx_corr]

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


ll_contributions <- function(demand_formula = NULL,
                             supply_formula = NULL,
                             data = NULL,
                             beta = NULL,
                             equal_sigmas = FALSE,                             na.action = na.exclude,
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


predict.diseq <- function(diseq_obj,
                          newdata = NULL,
                          type='prob_supply_constrained',
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
    orig_rownames <- diseq_obj$orig_rownames
  } else {
    list[X, coef_indices] = model.matrix.diseq(
      update(diseq_obj$demand_terms, NULL ~ .),
      update(diseq_obj$supply_terms, NULL ~ .),
      data = na.action(newdata[, union(all.vars(update(diseq_obj$demand_terms, NULL ~ .)),
                                       all.vars(update(diseq_obj$supply_terms, NULL ~ .))
                                       ),
                               with=FALSE]),
      corr = corr,
      equal_sigmas = equal_sigmas
      )
    orig_rownames <- rownames(newdata)
  }
  new_rownames <- rownames(X)

  idx_bd <- coef_indices[['beta_demand']]
  idx_bs <- coef_indices[['beta_supply']]
  idx_sd <- coef_indices[['sigma_demand']]
  idx_ss <- coef_indices[['sigma_supply']]
  idx_corr <- coef_indices[['sigma_corr']]
  
  beta_opt <- diseq_obj$raw_coefficients

  predicted_demand <- X[, idx_bd] %*% beta_opt[idx_bd]

  sigma1 <- exp(beta_opt[idx_sd])
  if (is.null(idx_ss)) {
    sigma2 <- sigma1
  } else {
    sigma2 <- exp(beta_opt[idx_ss])
  }
  if (is.null(idx_corr)) {
    rho <- 0
  } else {
    rho <- beta_opt[idx_corr]
  }

  if (type == 'prob_supply_constrained') {

    if (prob_without_positive_demand == TRUE) {

      sigma <- sqrt(sigma1^2 + sigma2^2 - 2 * corr * sigma1 * sigma2)
      prob_supply_constrained <- pnorm((predicted_demand-predicted_supply) / sigma)

    } else {

      mvnorm_fun <- if (exact == TRUE) {mvnorm_exact} else {mvnorm_approx}

      # A P(D > S metszet D > 0) = P(S - D < 0 metszet -D < 0) együttes eloszlása
      # levezetés papíron Peti mappájában (térkép!)
      mu_1 <- predicted_supply - predicted_demand
      mu_2 <- -predicted_demand
      sigma_1_bar <- sqrt(sigma1^2 + sigma2^2 - 2 * rho * sigma1 * sigma2)
      sigma_2_bar <- sigma2
      corr_coef <- sigma1 / sigma2 - rho

      prob_supply_constrained <- mvnorm_fun(mu_1, sigma_1_bar, mu_2, sigma_2_bar, corr_coef)

    }

    out <- prob_supply_constrained

  } else if (type == 'demand') {
    out <- predicted_demand
  } else if (type == 'supply') {
    out <- predicted_supply
  } else if (type == 'response') {
    out <- pmax(0, pmin(predicted_supply, predicted_demand))
  } else {
    stop('Unknown prediction type. Use one of: prob_supply_constrained, demand, supply, response.')
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


summary.diseq <- function(diseq_obj) {

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

  I_opt <- hessian(loglike, beta_opt)

  raw_std_err <- sqrt(diag(solve(I_opt)))
  raw_t_value <- beta_opt / raw_std_err
  raw_p_value <- pnorm(-abs(raw_t_value)) * 2

  coefficients <- diseq_obj$coefficients
  std_err <- raw_std_err
  std_err[idx_sd] <- std_err[idx_sd] * exp(beta_opt[idx_sd])  # Delta-method
  if (!is.null(idx_ss)) {
    std_err[idx_ss] <- std_err[idx_ss] * exp(beta_opt[idx_ss])  # Delta-method
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


export_to_csv <- function (diseq_summary_obj, file) {

  corr <- !is.null(diseq_summary_obj$settings$corr) && diseq_summary_obj$settings$corr == TRUE

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

