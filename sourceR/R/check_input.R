globalVariables(c("human", "value", "source_id")) 

## must be called after no_T, no_L and no_I are caclulated
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

check_mcmc_params <- function(mcmc_params, num) {
  if (!is.list(mcmc_params)) stop("mcmc_params must be a list.")
  ## save_lambda
  if (!"save_lambda" %in% names(mcmc_params)) {
    mcmc_params$save_lambda <- TRUE
  } else {
    if (!is.logical(mcmc_params$save_lambda)) stop("mcmc_params$save_lambda must be TRUE or FALSE.")
  }
  ## burn_in
  if (!"burn_in" %in% names(mcmc_params)) {
    mcmc_params$burn_in <- 0
  } else {
    if (!is.numeric(mcmc_params$burn_in) || !is.finite(mcmc_params$burn_in) || !is.wholenumber(mcmc_params$burn_in) || mcmc_params$burn_in < 0) 
      stop("mcmc_params$burn_in must be a whole positive number.")
  }
  ## thin
  if (!"thin" %in% names(mcmc_params)) {
    mcmc_params$thin <- 1
  } else {
    if (!is.numeric(mcmc_params$thin) || !is.finite(mcmc_params$thin) || !is.wholenumber(mcmc_params$thin)  || mcmc_params$thin <= 0) 
      stop("mcmc_params$thin must be a whole positive number.")
  }
  ## n_r
  if (!"n_r" %in% names(mcmc_params)) {
    mcmc_params$n_r <- ceiling(0.4 * num$no_I * num$no_T * num$no_J)
  } else {
    if (!is.numeric(mcmc_params$n_r) || !is.finite(mcmc_params$n_r) || !is.wholenumber(mcmc_params$n_r)  || mcmc_params$n_r <= 0) 
      stop("mcmc_params$n_r must be a whole positive number.")
  }
  return(mcmc_params)
}

check_params_fix <- function(params_fix) {
  if (!is.list(params_fix)) stop("params_fix must be a list.")
  ## r
  if (!"r" %in% names(params_fix)) {
    params_fix$r <- FALSE
  } else {
    if (!is.logical(params_fix$r)) stop("params_fix$r must be TRUE or FALSE.")
  }
  return(params_fix)
}

check_n_iter <- function(n_iter) {
  if (!is.numeric(n_iter) || !is.finite(n_iter) || !is.wholenumber(n_iter)  || n_iter <= 0) 
    stop("n_iter must be a positive whole number.")
  return(n_iter)
}

check_likelihood_dist <- function(likelihood_dist) {
  if (!is.character(likelihood_dist) || !(likelihood_dist == "pois" || likelihood_dist == "nbinom"))
    stop("likelihood_dist must be \"pois\" or \"nbinom\".")
  return(likelihood_dist)
}

check_and_process_data <- function(formula, type, time, location, data) {
  
  # !(class(time) == "formula") || !(class(location) == "formula")
  
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!is(formula, "formula")) stop("formula must be a formula.")
  if (missing (time)) {
    data$time <- rep(1, dim(data)[1])
    time_name <- "time"
  } else {
    if (!is(time, "formula")) stop("time must be a formula.")
    ## time column name
    time_name <- attributes(terms(time))$term.labels
  }
  if (missing (location)) {
    data$location <- rep("A", dim(data)[1])
    location_name <- "location"
  } else {
    if (!is(location, "formula")) stop("location must be a formula.")
    ## location column name
    location_name <- attributes(terms(location))$term.labels
    location_name <- gtools::mixedsort(location_name)
  }
  if (!is(type, "formula")) stop("type must be a formula.")
  
  matequal <- function(x, y) # used to make sure the source matrices are equal at each location as they can only vary over time
    is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)
  
  num <- list(no_I = NA, no_J = NA, no_L = NA, no_T = NA)
  
  ## human data name
  human_name <- rownames(attributes(terms(formula))$factors)[1]
  ## source names
  source_ids <- attributes(terms(formula))$term.labels
  source_ids <- gtools::mixedsort(source_ids)
  num$no_J <- length(source_ids)
  ## type column name
  type_name <- attributes(terms(type))$term.labels

  if (!time_name %in% names(data) || !location_name %in% names(data) || !human_name %in% names(data) || !all(source_ids %in% names(data))  || !type_name %in% names(data)) 
    stop("The data must be a data frame with columns labelled \"", human_name, "\", \"", paste(c(source_ids), collapse="\", \""), "\", \"", time_name, "\", \"", location_name, "\" and \"", type_name, "\". \n")
  
  data <- data[,c(human_name, source_ids, time_name, location_name, type_name)]
  names(data) <- c("human", source_ids, "time", "location", "type") # paste("source", 1 : num$no_J, sep = "")

  if (any(data[, c("human", source_ids)] < 0)) stop("All data must be positive.") # paste("source", 1 : num$no_J, sep = "")
  ## check all data is whole numbers
  is.whole <- function(x)
    is.numeric(x) && floor(x) == x 
  if (any(!apply(data[, c("human", source_ids)], c(1,2), is.whole))) stop("All data must be whole numbers.") # paste("source", 1 : num$no_J, sep = "")
  
  data$time <- as.factor(data$time)
  data$location <- as.factor(data$location)
  data$type <- as.factor(data$type)
  
  time_ids <- paste(gtools::mixedsort(unique(data$time)))
  location_ids <- paste(gtools::mixedsort(unique(data$location)))
  type_ids <- paste(gtools::mixedsort(unique(data$type)))
  
  num$no_L <- length(unique(data$location))
  num$no_T <- length(unique(data$time))
  num$no_I <- length(unique(data$type))
  
  if (dim(data)[1] != (num$no_L * num$no_T * num$no_I)) {
    stop("Some data is missing. The number of rows for each time and location combination should equal the number of types.")
  }
  
  data <- data[with(data, order(time, location, type)), ]
  time_ids <- unique(data$time)
  location_ids <- unique(data$location)
  
  source_data <- list()
  human_data <- list()
  r <- list()
  for (t in 1 : num$no_T) {
    human_data[[t]] <- list()
    for (l in 1 : num$no_L) {
      if (l != num$no_L) {
        ## check the source data is the same for each location within each time
        source_sub <- list()
        for (l1 in (l + 1) : num$no_L) {
          source_sub1 <- subset(data, subset = (time == time_ids[t] & location == location_ids[l]), select = c(-human, -time, -location))
          source_sub1 <- source_sub1[as.factor(type_ids), ]
          source_sub2 <- subset(data, subset = (time == time_ids[t] & location == location_ids[l1]), select = c(-human, -time, -location))
          source_sub2 <- source_sub2[as.factor(type_ids), ]
          if (! matequal(as.matrix(source_sub1), as.matrix(source_sub2))) 
            stop("The source data needs to be the same for each location within each time.")
        }
      }
      human_data[[t]][[l]] <- subset(data, subset = (time == time_ids[t] & location == location_ids[l]), select = c(human))
      rownames(human_data[[t]][[l]]) <- c(as.character(subset(data, subset = (time == time_ids[t] & location == location_ids[l]), select = c(type))$type))
      human_data[[t]][[l]]$human <- human_data[[t]][[l]][type_ids, "human"] # order rows and cols so that they are in alphabetical order
      
      names(human_data[[t]])[l] <- paste("location", location_ids[l], sep = "")
      names(human_data)[t] <- paste("time", time_ids[t], sep = "")
    }
    
    source_data[[t]] <- subset(data, subset = (time == time_ids[t] & location == location_ids[1]), select = c(-human, -time, -location, -type))
    rownames(source_data[[t]]) <- c(as.character(subset(data, subset = (time == time_ids[t] & location == location_ids[1]), select = c(type))$type))
    source_data[[t]] <- source_data[[t]][type_ids, source_ids] # order rows and cols so that they are in alphabetical order
    
    names(source_data)[t] <- paste("time", time_ids[t], sep = "")
    if (dim(source_data[[t]])[1] != num$no_I) stop("Some data is missing.") 
    else {
      source_data_no0 <- apply(source_data[[t]], c(1,2), function(ax) {if (ax == 0) 
        ax = 0.000001
      else 
        ax = ax
      return (ax)})
      
      r[[t]] <- apply(source_data_no0, 2, function(x) x / sum(x))
      names(r)[t] <- paste("time", time_ids[t], sep = "")
    }
  }
  return(list(r = r, source_data = source_data, human_data = human_data, num = num, data_names = list(source_ids = source_ids, type_ids = type_ids, time_ids = time_ids, location_ids = location_ids)))
}

check_priors <- function(priors, num, params_fix, likelihood_dist, formula, time, type, location, data_names) {
  if (!is.list(priors)) stop("priors must be a list.")

  ## check priors for source effects (accepts a single number)
  if (! "a" %in% names(priors)) stop("The parameters for the prior distributions for the source effects (a) must be specified.")
  else {
      if (missing(time) && missing(location)) {
        priors$a <- check_source_effects(data = priors$a, formula = formula, num = num, data_names = data_names)
      } else if (missing(time)) {
        priors$a <- check_source_effects(data = priors$a, formula = formula, location = location, num = num, data_names = data_names)
      } else if (missing(location)) {
        priors$a <- check_source_effects(data = priors$a, formula = formula, time = time, num = num, data_names = data_names)
      } else {
        priors$a <- check_source_effects(data = priors$a, formula = formula, time = time, location = location, num = num, data_names = data_names)
      }
  }
  
  ## check priors for r (accepts a single number or data frame containing a single number per time and source)
  if (! "r" %in% names(priors)) {
    if (params_fix$r) priors$r <- NA
    else stop("The parameters for the prior distribution for r must be specified.")
  } else {
    if (missing(time)) {
      if (is.data.frame(priors$r) == TRUE) {
        priors$r$time <- rep(1, dim(priors$r)[1])
      }
      priors$r <- check_r(priors$r, formula, type = type, num = num, data_names = data_names)
    } else {
      priors$r <- check_r(priors$r, formula, time, type, num, data_names)
    }
  }
  
  ## check priors for q (accepts a vector of length 2)
  if (! "theta" %in% names(priors)) {
    stop("The parameters for the base distribution for the type effects (theta) must be specified.")
  } else {
    if (! is.numeric(priors$theta) || length(priors$theta) != 2 || !all(is.finite(priors$theta)) || !all(priors$theta > 0)) 
      stop("The base distribution for the type effects (theta) must be a vector of length 2.")
  }
  
  ## check priors for d (accepts a vector of length 2)
  if (! "d" %in% names(priors)) {
    if (likelihood_dist == "pois") priors$d <- NA
    else stop("The parameters for the prior distribution for the dispersion parameter (d) \nmust be specified when using a Negative Binomial likelihood.")
  } else {
    if (! is.numeric(priors$d)) stop(" The prior for d must be a vector of length 2.")
    else if (length(priors$d) != 2 || !all(is.finite(priors$d)) || !all(priors$d > 0)) stop("The prior for d must be a vector of length 2.")
  }
  
  return(priors)
}

check_prev <- function(prev, formula, time, num, data_names) {
  
  if (is.data.frame(prev)) {
    if (missing(time)) {
      time_name = "time"
      prev$time <- rep(1, dim(prev)[1])
    } else {
      time_name <- attributes(terms(time))$term.labels
    }
    source_ids <- data_names$source_ids
    time_ids <- data_names$time_ids
    if (dim(prev)[1] != (num$no_T * num$no_J) || dim(prev)[2] != 3) {
      stop("The prevalence data must be a data frame of positive numbers \n
           with columns for the value, source_id and time. These must be named \"value\", \"source_id\" and \"", time_name, "\".\n")
    } else {
      if (is.numeric(prev$value) && is.finite(prev$value) && all(prev$value > 0)) {
        
        if (!time_name %in% names(prev) || !"value" %in% names(prev) || !"source_id" %in% names(prev)) 
          stop("The prevalence data must be a data frame of positive numbers \n
           with columns for the value, source_id and time. These must be named \"value\", \"source_id\" and \"", time_name, "\".\n")
        
        prev <- prev[,c("value", time_name, "source_id")]
        names(prev) <- c("value", "time", "source_id")
        
        source_ids_prev <- gtools::mixedsort(as.character(unique(prev$source_id)))
        if (!all.equal(source_ids_prev, source_ids)) stop("The factors in \"source_id\" must be the same as the source columns in the data data frame.")
        if (!all.equal(as.character(time_ids), as.character(unique(prev$time)))) stop("The factors in \"", time_name, "\" must be the same for the prev and data data frames.")
        
        prev <- prev[with(prev, order(time, source_id)), ]
        
        prev_list <- list()
        for (t in 1 : num$no_T) {
          prev_list[[t]] <- c(subset(prev, subset = (time == time_ids[t]), select = value)$value)
          if (any(prev_list[[t]] > 1) || any(prev_list[[t]] <= 0)) stop ("The prevalences must all be larger than 0 and smaller than or equal to 1.")
          names(prev_list)[t] <- paste("time", time_ids[t], sep = "")
        }
      } else stop("The prevalence data must be a data frame of positive numbers \n
           with columns for the value, source_id and time. These must be named \"value\", \"source_id\" and \"", time_name, "\".\n")
    }
  }
  return(prev_list)
}

check_initials <- function(initials, params_fix, likelihood_dist, priors, human_data, prev, num, r, formula, time, location, type, data_names) {
  if (!is.list(initials)) stop("initials must be a list.")
  
  ## initialise d
  if (likelihood_dist == "nbinom") {
    if (!"d" %in% names(initials)) initials$d <- rlnorm(1, priors$d[1], priors$d[2])
    else if (!is.numeric(initials$d) || length(initials$d) != 1 || initials$d <= 0) stop("The initial value for d must be a single positive number.")
  } else if (!"d" %in% names(initials)) initials$d <- 0
  
  ## initialise r
  if (!"r" %in% names(initials)) {
    initials$r <- r 
  } else {
      if (missing(time)) {
        if (is.data.frame(priors$r) == TRUE) {
          priors$r$time <- rep(1, dim(priors$r)[1])
        }
        initials$r <- check_r(initials$r, formula, type = type, num = num, data_names = data_names)
      } else {
        initials$r <- check_r(initials$r, formula, time, type, num, data_names)
      }
    for (t1 in 1 : num$no_T) {
      for (l1 in 1 : num$no_L) {
        for (sourcej in 1 : num$no_J) {
          sum_r = sum(initials$r[[t1]][[l1]][, sourcej])
          if (!all.equal(sum_r, 1)) { # all.equal() compares two objects using a numeric tolerance of .Machine$double.eps^0.5
            stop(paste("Error: the sum of the r initial values is ", sum_r, "for time ", t1, ", location ", l1, ", source ", colnames(initials$r)[sourcej], " rather than 1\n"))
        }
        }
      }
    }
  }
  
  ## initialise the source effects (a)
  if (!"a" %in% names(initials)) { 
    a_inits <- c(gtools::rdirichlet(1, rep(1, num$no_J * num$no_T * num$no_L)))
    initials$a <- list()
    k <- 0
    for (t in 1 : num$no_T) {
      initials$a[[t]] <- list()
      for (l in 1 : num$no_L) {
        k <- k + 1
        initials$a[[t]][[l]] <- a_inits[c((k * num$no_J - num$no_J + 1) : (k * num$no_J))]
      }
    }
  } else { # if initial values are provided
    if (missing(time) && missing(location)) {
      initials$a <- check_source_effects(data = initials$a, formula, num, data_names)
    } else if (missing(time)) {
      initials$a <- check_source_effects(data = initials$a, formula, location, num, data_names)
    } else if (missing(location)) {
      initials$a <- check_source_effects(data = initials$a, formula, time, num, data_names)
    } else {
      initials$a <- check_source_effects(data = initials$a, formula, time, location, num, data_names)
    }
    sum_a = sum(sapply(initials$a, function(x) sapply(lapply(x, '['), sum, simplify = T), simplify = T))
    if (!all.equal(sum_a, 1)) { # all.equal() compares two objects using a numeric tolerance of .Machine$double.eps^0.5
      stop(paste("Error: the sum of the source effects is ", sum_a, " rather than 1\n"))
    }
  }
  
  ## initialise alpha, q, theta and cluster (type effects related parameters)
  if (! "q" %in% names(initials)) { 
    ## use chinese restaurant process to simulate starting values for the q's 
    
    initialise_DP_vals <- function () {
      sum_const <- function(h) {
        sums <- 0
        for (t1 in 1 : num$no_T) {
          for (l1 in 1 : num$no_L) {
            ## calculate const for an individual h
            sums <- sums + sum(r[[t1]][h,, drop = F] %*% (initials$a[[t1]][[l1]] * prev[[t1]]))
          }
        }
        return (sums)
      }
      groups <- factor(rep(1, num$no_I))
      theta <- mean(unlist(human_data)) / mean(unlist(lapply(c(1 : num$no_I), function(x) sum_const(x))))
      qi <- rep(theta, num$no_I)
      
      ##checks
      stopifnot(length(theta) == length(table(groups)))
      stopifnot(length(groups) == num$no_I)
      stopifnot(length(qi) == num$no_I)
      
      return(list(cluster = groups, qi = qi, theta = theta))
    }
    
    res <- initialise_DP_vals()
    initials$q <- res$qi
    initials$cluster <- res$cluster
    initials$theta <- res$theta
  } else {
    if (is.data.frame(initials$q)) {
      if (!all(data_names$type_ids %in% names(initials$q$Type))) stop("The Type column must contain the same type names as those specified in the model formulation.")
      initials$q <- c(initials$q[data_names$type_ids, value])
      unique_q <- unique(initials$q)
      if (!is.numeric(unique_q) || length(unique_q) > num$no_I || any(unique_q <= 0)) {
        message("Warning: The number of unique values in the starting values of the type effects (q) must be positive and smaller than the maximum number of groups. \nThe starting values will be chosen for you.\n")
        res <- initialise_DP_vals()
        initials$q <- res$qi
        initials$cluster <- res$cluster
        initials$theta <- res$theta
      } else {
        if (likelihood_dist == "nbinom" && length(unique_q) > 1) 
          stop("Warning: All type effects must be in the same group and have the same starting value when using a negative binomial likelihood. \n")
        initials$theta <- unique_q     
        initials$cluster <- as.factor(as.numeric(as.factor(initials$q))) # as.numeric on a factor returns the level the factor is, not it's value (ie the 1st level of the factor...)
      }
    } else {
      message("Warning: initials$q must be a data frame with columns containing the initial values (named \"value\") and the corresponding type (named \"", data_names$type_ids, "\").\nThe starting values will be chosen for you.\n")
      res <- initialise_DP_vals()
      initials$q <- res$qi
      initials$cluster <- res$cluster
      initials$theta <- res$theta
    }
  }
  
  ## initialise lambda i
  initials$li <- calc_li(num$no_J, num$no_I, num$no_T, num$no_L, r = initials$r, a = initials$a, prev, q = initials$q)
  initials$lj <- calc_lj(num$no_J, num$no_I, num$no_T, num$no_L, r = initials$r, a = initials$a, prev, q = initials$q)
  
  return(initials)
}

check_alpha <- function(alpha_param, initials, likelihood_dist) {
  if (likelihood_dist == "pois") { # ie update type effects
    if (!is.numeric(alpha_param) || length(alpha_param) != 1 || alpha_param <= 0) stop("The initial value for alpha must be a single positive number.")
    initials$alpha <- alpha_param
    } else {
    initials$alpha <- NA
    }
  return(initials)
}

check_tuning <- function(tuning, num, params_fix, likelihood_dist, data_names) {
  if (!is.list(tuning)) stop("tuning must be a list.")
  
  adaptive_matrices <- list()
  
  ## no tuning for q, alpha, theta, cluster, pi and weights as they are Gibbs sampled
  
  ## No tuning for a: set up adaptive update
  n_update_a = list()
  n_accept_a = list()
  tuning_a = list()
  for (t in 1 : num$no_T) {
    n_update_a[[t]] = list()
    n_accept_a[[t]] = list()
    tuning_a[[t]] = list()
    for (l in 1 : num$no_L) {
      n_update_a[[t]][[l]] = rep(0, num$no_J)
      n_accept_a[[t]][[l]] = rep(0, num$no_J)
      tuning_a[[t]][[l]] = rep(0.01, num$no_J)
    }
  }
  
  adaptive_matrices$a <- list()
  adaptive_matrices$a$n_update_a <- n_update_a
  adaptive_matrices$a$n_accept_a <- n_accept_a
  adaptive_matrices$a$tuning_a <- tuning_a
  
  ## No tuning for r: set up adaptive update
  # set up matrices for adaptive updates
  # Every time a cell of r is proposed to be updated, add 1 to n_update_r[i, j]
  # If accepted add 1 to n_accept_r[i, j]
  # Every 50 or so iterations check the acceptance rate (n_accept_r[i, j] / n_update_r[i, j])
  # If it is below 0.45 decrease tuning parameter by delta(n) = min(0.01, n ^ -1/2), else increase by the same
  # Diminishing adaptation as n increases
  # Bounded by 0.01 (so it can't adapt too much)
  n_update_r <- list()
  n_accept_r <- list()
  tuning_r <- list()
  for (t in 1 : num$no_T) {
    n_update_r[[t]] <- matrix(0, ncol = num$no_J, nrow = num$no_I)
    n_accept_r[[t]] <- matrix(0, ncol = num$no_J, nrow = num$no_I)
    tuning_r[[t]] <- matrix(0.01, ncol = num$no_J, nrow = num$no_I)
  }
  adaptive_matrices$r <- list()
  adaptive_matrices$r$n_update_r <- n_update_r
  adaptive_matrices$r$n_accept_r <- n_accept_r
  adaptive_matrices$r$tuning_r <- tuning_r
  
  ## No tuning for d: set up adaptive update
  ## Same algorithm as for r
  if (likelihood_dist == "nbinom") {
    #   if (! "d" %in% names(tuning)) {
    #     fixed_d <- 0.01
    #   } else {
    #     if (!is.numeric(fixed_d) || length(fixed_d) != 1 || !is.finite(fixed_d) || fixed_d <= 0) stop("The fixed tuning value for d must be a numeric vector of length 1.")
    #     else fixed_d <- tuning$d
    #   }
    adaptive_matrices$d <- list()
    adaptive_matrices$d$n_update_d <- 0
    adaptive_matrices$d$n_accept_d <- 0
    adaptive_matrices$d$tuning_d <- 0.2
  } else {
    # fixed_d <- NA
    adaptive_matrices$d <- NA
  }
  
  # aosp gives the adaptive optimal scaling parameter for block adaptive updates
  # 2.38 is optimal for Gaussian updates
  return(adaptive_matrices = adaptive_matrices)
}

create_posterior <- function(num, likelihood_dist, params_fix, mcmc_params, data_names, initials, r, n_iter) {
  posterior <- list()
  
  ## calculate the number of iterations that will be saved
  save_i <- seq(mcmc_params$burn_in + 1, n_iter, by = mcmc_params$thin)
  n_iter_thinned <- length(save_i)
  n_iter <- max(save_i)
  iter_val <- 1
  
  params_cur <- initials
  params_can <- initials
  
  if (params_fix$r == FALSE) posterior$r <- list()
  if (likelihood_dist == "pois") {
    posterior$q <- matrix(NA, nrow = n_iter_thinned, ncol = num$no_I) 
    colnames(posterior$q) <- data_names$type_ids #paste("type", data_names$type_ids, sep = "")
    posterior$cluster <- matrix(NA, nrow = n_iter_thinned, ncol = num$no_I)
    posterior$theta <- matrix(NA, nrow = n_iter_thinned, ncol = num$no_I)
  } else {
    posterior$q <- matrix(NA, nrow = n_iter_thinned, ncol = 1) 
    colnames(posterior$q) <- "overall_mean"
  }

  if (likelihood_dist == "nbinom") {
    posterior$d <- rep(NA, n_iter)
    posterior$d[1] <- initials$d
  } 
  if (mcmc_params$save_lambda == T) {
    posterior$li <- list()
    posterior$lj <- list()
  }
  
  posterior$a <- list()
  for (t in 1 : num$no_T) {
    posterior$a[[t]] <- list()
    if (mcmc_params$save_lambda == T) {
      posterior$li[[t]] <- list()
      posterior$lj[[t]] <- list()
    }
    if (params_fix$r == F) {
      posterior$r[[t]] <- array(NA, dim = c(num$no_I, num$no_J, n_iter_thinned))
      dimnames(posterior$r[[t]]) <- list(data_names$type_ids, data_names$source_ids, 1 : n_iter_thinned)
      posterior$r[[t]][, , 1] <- params_cur$r[[t]]
    }
    for (l in 1 : num$no_L) {
      posterior$a[[t]][[l]] <- matrix(NA, nrow = n_iter_thinned, ncol = num$no_J)
      colnames(posterior$a[[t]][[l]]) <- data_names$source_ids
      posterior$a[[t]][[l]][1, ] <- params_cur$a[[t]][[l]]
      if (mcmc_params$save_lambda == T) {
        posterior$li[[t]][[l]] <- matrix(NA, nrow = n_iter_thinned, ncol = num$no_I)
        colnames(posterior$li[[t]][[l]]) <- data_names$type_ids #paste("type", data_names$type_ids, sep = "")
        posterior$lj[[t]][[l]] <- matrix(NA, nrow = n_iter_thinned, ncol = num$no_J)
        colnames(posterior$lj[[t]][[l]]) <- data_names$source_ids
      }
    }
  }

  if (params_fix$r == FALSE) {
    names(posterior$r) <- paste("time", data_names$time_ids, sep = "")
  }
  if (mcmc_params$save_lambda == TRUE) {
    names(posterior$li) <- paste("time", data_names$time_ids, sep = "")
    names(posterior$lj) <- paste("time", data_names$time_ids, sep = "")
  }
  names(posterior$a) <- paste("time", data_names$time_ids, sep = "")
  
  for (t in 1 : num$no_T) {
    names(posterior$a[[t]]) <- paste("location", data_names$location_ids, sep = "")
    if (mcmc_params$save_lambda == TRUE) {
      names(posterior$li[[t]]) <- paste("location", data_names$location_ids, sep = "")
      names(posterior$lj[[t]]) <- paste("location", data_names$location_ids, sep = "")
    }
  }

  return(list(posterior = posterior, params_cur = params_cur, params_can = params_can))
}

create_acceptance_rate <- function(num, params_fix, likelihood_dist) {
  acceptance <- list()
  acceptance_rate <- list()
  n_accept <- list()
  n_update <- list()
  
  n_accept$a <- list()
  n_accept$r <- list()
  
  n_update$a <- list()
  n_update$r <- list()
  
  acceptance$a <- list()
  acceptance$r <- list()
  
  acceptance_rate$a <- list()
  if (params_fix$r == FALSE) acceptance_rate$r <- list()
  
  for (t in 1 : num$no_T) {
    n_accept$a[[t]] <- list()
    n_accept$r[[t]] <- 0
    
    n_update$a[[t]] <- list()
    n_update$r[[t]] <- 0
    
    acceptance$a[[t]] <- list()
    acceptance$r[[t]] <- 0
    
    acceptance_rate$a[[t]] <- list()
    acceptance_rate$r[[t]] <- 0
    for (l in 1 : num$no_L) {
      n_accept$a[[t]][[l]] <- numeric(num$no_J)
      
      n_update$a[[t]][[l]] <- numeric(num$no_J)
      
      acceptance$a[[t]][[l]] <- numeric(num$no_J)
      
      acceptance_rate$a[[t]][[l]] <- numeric(num$no_J)
    }
  }
  acceptance$d <- 0
  
  if (likelihood_dist == "nbinom") acceptance_rate$d <- 0
  
  return(list(acceptance = acceptance, acceptance_rate = acceptance_rate, n_accept = n_accept, n_update = n_update))
}

## used for checking the priors and initials
check_source_effects <- function(data, formula, time, location, num, data_names) { 
  ## checks/creates initials and prior matrices for the source effects

  source_ids <- data_names$source_ids
  time_ids <- data_names$time_ids
  location_ids <- data_names$location_ids
  
  if (is.data.frame(data)) {
    if (missing(time)) {
      time_name <- "time"
      data$time <- rep(1, dim(data)[1])
    } else {
      time_name <- attributes(terms(time))$term.labels
    }
    
    if (missing(location)) {
      location_name <- "location"
      data$location <- rep("A", dim(data)[1])
    } else {
      location_name <- attributes(terms(location))$term.labels
    }
    if (dim(data)[1] != (num$no_J * num$no_T * num$no_L) || dim(data)[2] != 4) {
      stop("The priors, initials and tuning values for the source effects must be a single positive number or a data frame of single positive numbers \n
              with columns for the value, time, location and source ids. \n
              These must be named \"value\", \"", time_name, "\", \"", location_name, "\", and \"source_id\".\n
              The factors in the \"source_id\" column must be the same as the names of the source columns in the data data frame.")
    } else {
      if (!time_name %in% names(data) || !location_name %in% names(data) || !"source_id" %in% names(data) || !"value" %in% names(data) ) 
        stop("The priors, initials and tuning values for the source effects must be a single positive number or a data frame of single positive numbers \n
              with columns for the value, time, location and source ids. \n
              These must be named \"value\", \"", time_name, "\", \"", location_name, "\", and \"source_id\".\n
              The factors in the \"source_id\" column must be the same as the names of the source columns in the data data frame.")
      
      data <- data[,c("value", "source_id", time_name, location_name)]
      names(data) <- c("value", "source_id", "time", "location")
      data <- data[with(data, order(time, location, source_id)), ]
      
      if (is.numeric(data$value) != TRUE || !all(is.finite(data$value)) || !all(data$value > 0)) {
        stop("The initial values for the source effects must be a single positive number or a data frame of single positive numbers \n
              with columns for the value, time, location and source ids. \n
              These must be named \"value\", \"", time_name, "\", \"", location_name, "\", and \"source_id\".\n
              The factors in the \"source_id\" column must be the same as the names of the source columns in the data data frame.")
      }
      
      if (!all.equal(gtools::mixedsort(as.character(unique(data$source_id))), source_ids)) stop("The factors in \"source_id\" must be the same as the source columns in the data data frame.")
      if (!all.equal(as.character(time_ids), gtools::mixedsort(as.character(unique(data$time))))) stop("The factors in \"", time_name, "\" must be the same for the prior, tuning, initial and data data frames.")
      if (!all.equal(as.character(location_ids), gtools::mixedsort(as.character(unique(data$location))))) stop("The factors in \"", location_name, "\" must be the same for the prior, tuning, initial and data data frames.")
      
      data_list <- list()
      for (t in 1 : num$no_T) {
        data_list[[t]] <- list()
        for (l in 1 : num$no_L) {
          data_list[[t]][[l]] <- c(subset(data, subset = (time == time_ids[t] & location == location_ids[l]), select = value))$value
          data_list[[t]][[l]] <- data_list[[t]][[l]][order(as.character(subset(data, subset = (time == time_ids[t] & location == location_ids[l]), select = source_id)$source_id))]
          names(data_list[[t]])[l] <- paste("location", location_ids[l], sep = "")
        }
        names(data_list)[t] <- paste("time", time_ids[t], sep = "")
      }
    }
  } else {
    if (length(data) != 1 || is.numeric(data) != TRUE || !all(is.finite(data)) || !all(data > 0)) {
      stop("The initial values for the source effects must be a single positive number or a data frame of single positive numbers \n
              with columns for the value, time, location and source ids. \n
              These must be named \"value\", \"", time_name, "\", \"", location_name, "\", and \"source_id\".\n
              The factors in the \"source_id\" column must be the same as the names of the source columns in the data data frame.")
    } else {
      # create the nested list of priors from the vector (same priors for each time and location)
        data_list <- list()
        for (t in 1 : num$no_T) {
          data_list[[t]] <- list()
          for (l in 1 : num$no_L) {
            data_list[[t]][[l]] <- rep(data, num$no_J)
            names(data_list[[t]])[l] <- paste("location", location_ids[l], sep = "")
          }
          names(data_list)[t] <- paste("time", time_ids[t], sep = "")
        }
    }
  }
  
  return(data_list)
}

## only used for checking the prior for r
check_r <- function(data, formula, time, type, num, data_names) { # data = priors$r
 
  if (is.data.frame(data)) {
    if (missing(time)) {
      time_name <- "time"
      data$time <- rep(1, dim(data)[1])
    } else {
      time_name <- attributes(terms(time))$term.labels
    }
    type_name <- attributes(terms(type))$term.labels
    
    if (dim(data)[1] != num$no_T * num$no_I * num$no_J || dim(data)[2] != 4) {
      stop("The tuning and prior values for the r must be a single positive number or a data frame of single positive numbers \n
           with columns for the value, time, type and source ids. \n
              These must be named \"value\", \"", time_name, "\", \"", type_name, "\", and \"source_id\".\n
              The factors in the \"source_id\" column must be the same as the names of the source columns in the data data frame.")
    } else {
      if (!time_name %in% names(data) || !type_name %in% names(data) || !"source_id" %in% names(data) || !"value" %in% names(data) ) 
        stop("The tuning and prior values for the r must be a single positive number or a data frame of single positive numbers \n
           with columns for the value, time, type and source ids. \n
              These must be named \"value\", \"", time_name, "\", \"", type_name, "\", and \"source_id\".\n
              The factors in the \"source_id\" column must be the same as the names of the source columns in the data data frame.")
      
      data <- data[,c("value", "source_id", time_name, type_name)]
      names(data) <- c("value", "source_id", "time", "type")
      data <- data[with(data, order(time, type, source_id)), ]
      time_ids <- data_names$time_ids
      type_ids <- data_names$type_ids
      source_ids <- data_names$source_ids
   
      if (is.numeric(data$value) != TRUE || !all(is.finite(data$value)) || !all(data$value > 0)) {
        stop("The tuning and prior values for the r must be a single positive number or a data frame of single positive numbers \n
           with columns for the value, time, type and source ids. \n
              These must be named \"value\", \"", time_name, "\", \"", type_name, "\", and \"source_id\".\n
              The factors in the \"source_id\" column must be the same as the names of the source columns in the data data frame.")
      } else {
        
        if (!all.equal(source_ids, gtools::mixedsort(as.character(unique(data$source_id))))) stop("The factors in \"source_id\" must be the same as the source columns in the data dataframe.")
        if (!all.equal(gtools::mixedsort(as.character(time_ids)), gtools::mixedsort(as.character(unique(data$time))))) stop("The factors in \"", time_name, "\" must be the same for the prior, initial and data dataframes.")
        if (!all.equal(gtools::mixedsort(as.character(type_ids)), gtools::mixedsort(as.character(unique(data$type))))) stop("The factors in \"", type_name, "\" must be the same for the prior, initial and data dataframes.")
        
        data_list <- list()
        for (t in 1 : num$no_T) {
          data_list[[t]] <- reshape(data[which(data$time == time_ids[t]), c("value", "source_id", "type")], timevar = "source_id", idvar = "type", direction = "wide")
          data_list[[t]] <- data_list[[t]][as.factor(gtools::mixedsort(as.character(data_list[[t]]$type))), ]
          rownames(data_list[[t]]) <- data_list[[1]]$type
          data_list[[t]] <- data_list[[t]][,-1]
          colnames(data_list[[t]]) <- source_ids
        }
      }
    }
  } else if (is.numeric(data) && length(data) == 1 && is.finite(data) && all(data > 0)) {
    data_list <- list()
    for (t in 1 : num$no_T) {
      data_list[[t]] <- matrix(rep(data, num$no_I * num$no_J), ncol = num$no_J)
    }
  } else stop("The priors for the r must be a single positive number or a data frame of single positive numbers \n
              with columns for the value and time. These must be named \"value\" and \"", time_name, "\".\n")
  
  return(data_list)
}
