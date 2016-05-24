#' @export
saBayes <- function(formula, time, location, type, data, priors, alpha_conc, prev, likelihood_dist, n_iter, mcmc_params, initials, params_fix) {
  ### data is long form
  ### formular must be of the form: y ~ x1 + ... + xJ 
  ##   (where y is the name of the vector of human cases, and x1,...,xJ are the names of the vectors of source counts)
  ### time is a vector of times. This is used to seperate the r matrix and the human cases into different times and use different (independent) 
  ##   source effects for each time
  ### location is a vector of locations. This is used to seperate the human cases into seperate locations and use a different 
  ##   (independent) source effect for each location.
  ### The source counts do not vary for different locations, hence the source cases must be repeated for each location.
  ### There must be the same number of types (ie rows) within each time and location combination. 
  
  ### mcmc_params is a list of parameters to control the speed and memory consumption of the algorithm
  ##   It contains a list with items named: save_lambda, gmax, n_r, burnin and thin.
  ##   Defaults: save_lambda (TRUE), burn_in (0), thin (1), and n_r (0.4 * number of locations * number of times)
  
  ### values for tuning and initial values are calculated automatically if not provided
  ### priors and initials are given in long form
  
  ## C++ code to calculate lambda i #################################################################### 
  calc_li <- function(no_J, no_I, no_T, no_L, r, a, prev, q, g) {
    .Call('sourceR_calc_li', PACKAGE = 'sourceR', no_J, no_I, no_T, no_L, r, a, prev, q, g)
  }
  
  ## C++ code to calculate lambda j ####################################################################
  calc_lj <- function(no_J, no_I, no_T, no_L, r, a, prev, q, g) {
    .Call('sourceR_calc_lj', PACKAGE = 'sourceR', no_J, no_I, no_T, no_L, r, a, prev, q, g)
  }
 
  ## check likelihood_dist
  if (missing(likelihood_dist)) stop("likelihood_dist is missing.")
  else likelihood_dist <- check_likelihood_dist(likelihood_dist)
  
  ## check n_iter
  if (missing(n_iter)) stop("n_iter is missing.")
  else n_iter <- check_n_iter(n_iter)
  
  ## Create params_fix
  # a1 should not be fixed when recaling the type effects to have mean 1 without using an overall mean
  # an overall mean is not implemented, so a1 is always fixed at present.
  if (missing(params_fix)) {
    if (likelihood_dist == "pois") params_fix <- list(type_effects = FALSE, a1 = FALSE, r = FALSE, alpha = FALSE)
    else params_fix <- list(type_effects = TRUE, a1 = FALSE, r = FALSE, alpha = TRUE)
  }
  else params_fix <- check_params_fix(params_fix)
  
  ## check data and formula
  if (missing(data)) stop("data is missing.") 
  else if (!(is.data.frame(data) || is.matrix(data))) stop("data must be a data frame or matrix.")
  data <- as.data.frame(data)
  
  if (missing(formula) || missing(type)) stop("The formula and/ or type is missing.")
  else if (!(class(formula) == "formula") || !(class(type) == "formula")) stop("formula, type, time and location must be valid formula objects.")
  else {
    if (missing(time) && missing(location)) {
      res <- check_and_process_data(formula = formula, type = type, data = data)
    } else if (missing(time)) {
      res <- check_and_process_data(formula = formula, type = type, location = location, data = data)
    } else if (missing(location)) {
      res <- check_and_process_data(formula = formula, type = type, time, data = data)
    } else {
      res <- check_and_process_data(formula = formula, type = type, time, location = location, data = data)
    }
  }
 
  r <- res$r
  source_data <- res$source_data
  human_data <- res$human_data
  num <- res$num
  data_names <- res$data_names
  
  ### check priors
  if (missing(priors)) stop("priors is missing.")
  else {
    if (missing(time) && missing(location)) {
      priors <- check_priors(priors, num, params_fix, likelihood_dist, formula, type = type, data_names = data_names)
    } else if (missing(time)) {
      priors <- check_priors(priors, num, params_fix, likelihood_dist, formula, type = type, location = location, data_names = data_names)
    } else if (missing(location)) {
      priors <- check_priors(priors, num, params_fix, likelihood_dist, formula, time = time, type = type, data_names = data_names)
    } else {
      priors <- check_priors(priors, num, params_fix, likelihood_dist, formula, time = time, type = type, location = location, data_names = data_names)
    }
  }
  
  ## check prev
  if (missing(prev)) { # if no prevalences supplied
    warning("Warning: Prevalence data for each source has not been supplied. \nA prevalence of 1 will be used for each source.\n")
    prev <- list()
    for (t in 1 : num$no_T) {
      prev[[t]] <- c(rep(1, num$no_J)) 
    }
  } else { # if a vector or lists of vectors supplied
    if (missing(time)) {
      prev <- check_prev(prev, formula, num = num, data_names = data_names)
    } else {
      prev <- check_prev(prev, formula, time, num, data_names)
    }
  }

  ## check mcmc params
  if (missing(mcmc_params)) {
    mcmc_params <- list()
  }
  mcmc_params <- check_mcmc_params(mcmc_params, num)

  ## check and create initials
  if (missing(initials)) {
    initials <- list()
  } 
  if (missing(time) && missing(location)) {
    initials <- check_initials(initials, params_fix, likelihood_dist, priors, human_data, prev, num, r, formula, type, data_names)
  } else if (missing(time)) {
    initials <- check_initials(initials, params_fix, likelihood_dist, priors, human_data, prev, num, r, formula, location, type, data_names)
  } else if (missing(location)) {
    initials <- check_initials(initials, params_fix, likelihood_dist, priors, human_data, prev, num, r, formula, time, type, data_names)
  } else {
    initials <- check_initials(initials, params_fix, likelihood_dist, priors, human_data, prev, num, r, formula, time, location, type, data_names)
  }

  ## check alpha
  if (missing(alpha_conc)) stop("alpha_conc is missing.")
  initials <- check_alpha(alpha_conc, initials, likelihood_dist)
  
  ## check and create tuning
  tuning <- list()
  tuning <- check_tuning(tuning, num, params_fix, likelihood_dist, data_names)
  adaptive_matrices <- tuning

  res <- create_posterior(num, likelihood_dist, params_fix, mcmc_params, data_names, initials, r, n_iter)
  posterior = res$posterior
  params_cur = res$params_cur
  params_can = res$params_can

  res <- doMCMC(data = list(source_data = source_data, human_data = human_data, r = r, prev = prev), num, priors, initials, likelihood_dist, 
         n_iter, posterior, params_cur, adaptive_matrices, mcmc_params, params_fix, data_names)
  return(res)
}





