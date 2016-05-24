disclapmix <- function(x, clusters, init_y = NULL, iterations = 100L, eps = 0.001, verbose = 0L, 
  glm_method = "internal_coef", glm_control_maxit = 50L, glm_control_eps = 1e-6, init_y_method = "pam", ...) {
  
  dots <- list(...)
  
  if ("centers" %in% names(dots)) {
    stop("centers argument has been deprecated, please refer to the manual. You can use your prefered parallelisation strategy.")
  }
  
  if ("use.parallel" %in% names(dots)) {
    stop("use.parallel argument has been deprecated, please refer to the manual. You can use your prefered parallelisation strategy.")
  }
  
  if ("calculate.logLs" %in% names(dots)) {
    stop("calculate.logLs argument has been deprecated, please refer to the manual. You can use verbose = 2L to get the log likelihoods during the iterations.")
  }
  
  if ("plots.prefix" %in% names(dots)) {
    stop("plots.prefix argument has been deprecated, please refer to the manual.")
  }
    
  ##############################################################################
  if (is.null(clusters) || length(clusters) != 1L) {
    stop("clusters must be a number, not a vector of numbers")
  }
      
  if (!is.integer(clusters)) {
    stop("clusters must be an integer, e.g. 1L, 5L or similar (note the required L suffix)")
  }
      
  if (clusters < 1L) {
    stop("clusters must be at least 1L")
  }
  
  ##
  
  if (is.null(verbose) || !is.integer(verbose) || length(verbose) != 1 || verbose < 0L | verbose > 2L) {
    stop("verbose must be an integer between 0L and 2L (inclusive, note the required L suffix)")
  }
  
  if (is.null(iterations) || length(iterations) != 1 || !is.integer(iterations) || iterations < 0L) {
    stop("iterations must be an integer strictly greater than 0L (note the required L suffix)")
  }
  
  if (is.null(glm_method) | !is.character(glm_method) | length(glm_method) != 1L | (glm_method != "glm.fit" & glm_method != "internal_dev" & glm_method != "internal_coef")) {
    stop("For now, only 'glm.fit', 'internal_dev' and 'internal_coef' are supported for glm_method.")
  }
  
  ##
  
  if (is.null(glm_control_maxit) || length(glm_control_maxit) != 1L) {
    stop("glm_control_maxit must be a number, not a vector of numbers")
  }
      
  if (!is.integer(glm_control_maxit)) {
    stop("glm_control_maxit must be an integer, e.g. 25L or similar (note the required L suffix)")
  }
      
  if (glm_control_maxit < 1L) {
    stop("glm_control_maxit must be at least 1L")
  }
  
  ##
  
  if (is.null(glm_control_eps) || length(glm_control_eps) != 1L) {
    stop("glm_control_eps must be a number, not a vector of numbers")
  }
      
  if (!is.numeric(glm_control_eps)) {
    stop("glm_control_eps must be a number, e.g. 1e-4 or similar")
  }
      
  if (glm_control_eps <= 0) {
    stop("glm_control_eps must be greater than 0")
  }
  
  ##
  
  if (is.null(init_y) && (
        is.null(init_y_method) || !is.character(init_y_method) || length(init_y_method) != 1L || (init_y_method != "pam" && init_y_method != "clara")
     )) {
    stop("The specified init_y_method is not valid, please refer to the documentation.")
  }
  
  if (!is.null(init_y) && !is.null(init_y_method)) {
    warning(paste("A init_y_method specified, '", init_y_method, "', will be ignored as init_y is supplied.", sep = ""))
    init_y_method <- NULL
  }

  check_x(x)
  
  #if (glm_method == "glm.fit" && clusters == 1L && ncol(x) == 1L) {
  #  warning("Only one cluster and one locus, using the internal_dev method.")
  #  glm_method <- "internal_dev"
  #}
  
  if (verbose >= 1L) {
    cat(as.character(Sys.time()), ": Starting estimation for ", clusters, " clusters.\n", sep = "")
  }
  
  y <- NULL
  
  # y
  if (!is.null(init_y)) {
    check_y(init_y, clusters = clusters, loci = ncol(x))
    
    y <- init_y

    if (verbose >= 1L) {
      cat(as.character(Sys.time()), ": Using supplied y as the initial central haplotypes.\n", sep = "")
    }
  } else {
    if (verbose >= 1L) {
      cat(as.character(Sys.time()), ": Estimating initial central haplotypes, y, using ", init_y_method, ".\n", sep = "")
    }
    
    y <- create_initial_y(x, clusters = clusters, init_y_method = init_y_method)
  
    if (verbose >= 1L) {
      cat(as.character(Sys.time()), ": Initial central haplotypes, y, estimated.\n", sep = "")
    }
  }

  if (verbose >= 1L) {
    if (glm_method == "internal_dev" || glm_method == "internal_coef") {
      cat(as.character(Sys.time()), ": No need to generate model matrix when using ", glm_method, ".\n", sep = "")
    } else {
      cat(as.character(Sys.time()), ": Starting to generate model matrix.\n", sep = "")
    }
  }
  
  model.matrix.function <- model.matrix
 
  #if (glm_method == "OLD_internal" || glm_method == "internal") {
  #  model.matrix.function <- sparse.model.matrix
  #}
  
  model_matrix <- NULL
  
  if (!(glm_method == "internal_dev" || glm_method == "internal_coef")) {
    model_matrix <- create_model_matrix(x, clusters, model.matrix.function)
    
    if (verbose >= 1L) {
      cat(as.character(Sys.time()), ": Model matrix generated.\n", sep = "")
    }
  }
  
  if (verbose >= 1L) {
    cat(as.character(Sys.time()), ": Starting to generate initial response vector.\n", sep = "")
  }
  
  response_vector <- create_response_vector(x, y)
  
  if (verbose >= 1L) {
    cat(as.character(Sys.time()), ": Initial response vector done.\n", sep = "")
  }
  
  vic_matrix <- matrix(1/clusters, nrow = nrow(x), ncol = nrow(y))
  #vic_matrix <- matrix(runif(nrow(x)*nrow(y), 0.2, 0.8), nrow = nrow(x), ncol = nrow(y))
  #vic_matrix <- vic_matrix / rowSums(vic_matrix)
  weight_vector <- rcpp_create_new_weight_vector(vic_matrix, ncol(y))
  tau_vector <- rep(1/clusters, clusters)

  ##############################################################################
  if (verbose >= 1L) {
    cat(as.character(Sys.time()), ": Model matrix and initial vectors has been generated.\n", sep = "")
  }
  ##############################################################################
  
  disclap_parameters <- NULL
  
  fit <- NULL

  iterations_total <- 0L
  
  logL_full_iterations <- NULL
  logL_marginal_iterations <- NULL
  AIC_full_iterations <- NULL
  AIC_marginal_iterations <- NULL
  AICc_full_iterations <- NULL
  AICc_marginal_iterations <- NULL
  BIC_full_iterations <- NULL
  BIC_marginal_iterations <- NULL
  
  v_gain_iterations <- NULL
  tau_iterations <- tau_vector
  changed_center <- NULL
  centers_iterations <- list()

  converged <- FALSE

  if (verbose >= 1L) {
    cat(as.character(Sys.time()), ": Starting the EM algorithm using ", glm_method, " IRLS method.\n", sep = "")
  }
  
  ##############################################################################
  if (clusters == 1L) {
    # Only one cluster, only one iteration in the EM algorithm is 
    # necessary (the weights are all 1 and no reestimation is needed)
    iterations <- 1L
    converged <- TRUE
  }
  
  repeat {   
    last_v_max <- 0
    last_logL_full <- -Inf
    last_logL_marginal <- -Inf 

    ############################################################################
    # EM
    ############################################################################
    for (iter in 1L:iterations) {
      iterations_total <- iterations_total + 1L
      
      #################################
      # M-step 
      #################################
      fit <- NULL
      beta_start <- NULL
      
      if (iter == 1L && nrow(y) > 1L) {
        beta_start <- c(rep(-1, nrow(y)), rep(0.5, ncol(x)-1))
        #print(beta_start)
        #print(head(model_matrix))          
      }
      
      disclap_parameters <- NULL
      
      #print(head(vic_matrix))
 
      if (glm_method == "glm.fit") {
        fit <- glm.fit(y = response_vector, x = model_matrix, 
          intercept = FALSE, weights = weight_vector, family = DiscreteLaplace(),
          start = beta_start,
          control = glm.control(trace = (verbose >= 2L), epsilon = glm_control_eps, maxit = glm_control_maxit))
        disclap_parameters <- convert_coef_to_disclap_parameters(fit$coefficients, clusters)
      } else if (glm_method == "internal_dev") {
        fit <- INTERNAL_glmfit(loci = ncol(x), clusters = clusters, individuals = nrow(x), 
          response_vector = response_vector, apriori_probs = tau_vector, weight_vector = weight_vector, vmat = vic_matrix, 
          verbose = (verbose >= 2L), stop_by_deviance = TRUE, epsilon = glm_control_eps, maxit = glm_control_maxit) 
        disclap_parameters <- convert_coef_to_disclap_parameters_internal(fit$coefficients, clusters)
        colnames(disclap_parameters) <- colnames(x)
      } else if (glm_method == "internal_coef") {
        fit <- INTERNAL_glmfit(loci = ncol(x), clusters = clusters, individuals = nrow(x), 
          response_vector = response_vector, apriori_probs = tau_vector, weight_vector = weight_vector, vmat = vic_matrix, 
          verbose = (verbose >= 2L), stop_by_deviance = FALSE, epsilon = glm_control_eps, maxit = glm_control_maxit) 
        disclap_parameters <- convert_coef_to_disclap_parameters_internal(fit$coefficients, clusters)
        colnames(disclap_parameters) <- colnames(x)
      } else {
        stop("Unsupported glm_method chosen")
      }
      
      rownames(disclap_parameters) <- paste("cluster", 1L:clusters, sep = "")    
      
      #print(fit$coefficients)
      #print(disclap_parameters)

      if (fit$converged == FALSE) {
        msg <- paste(glm_method, " IRLS did not converge in iteration ", iterations_total, sep = "")
        warning(msg)
      }
      
      #################################
      # E-step ########################
      #################################
      wic <- rcpp_calculate_wic(x, y, disclap_parameters, tau_vector)
      
      new_vic_matrix <- rcpp_calculate_vic(wic)
      #print(head(new_vic_matrix))
      new_tau_vector <- apply(new_vic_matrix, 2, sum) / nrow(x)
      new_weight_vector <- rcpp_create_new_weight_vector(new_vic_matrix, ncol(y))
      tau_iterations <- rbind(tau_iterations, new_tau_vector)
      
      #print(disclap_parameters)
      #print(new_tau_vector)
      #################################
      # v gain ########################
      #################################
      v_max_diff <- max(abs(new_vic_matrix - vic_matrix))
      v_gain <- v_max_diff / last_v_max
      v_gain_iterations <- c(v_gain_iterations, v_gain)
      
      #################################
      # Updating ######################
      #################################
      vic_matrix <- new_vic_matrix
      tau_vector <- new_tau_vector
      weight_vector <- new_weight_vector
      
      if (verbose >= 2L) {
        logL_full <- get_loglikelihood_full(fit, clusters, ncol(x), response_vector, weight_vector, tau_norm = tau_vector^(1/ncol(x)))
        logL_marginal <- get_loglikelihood_marginal(x, y, disclap_parameters, tau_vector)
        AIC_full <- get_AIC(logL_full, individuals = nrow(x), clusters = clusters, loci = ncol(x))
        AIC_marginal <- get_AIC(logL_marginal, individuals = nrow(x), clusters = clusters, loci = ncol(x))
        AICc_full <- get_AICc(logL_full, individuals = nrow(x), clusters = clusters, loci = ncol(x))
        AICc_marginal <- get_AICc(logL_marginal, individuals = nrow(x), clusters = clusters, loci = ncol(x))
        BIC_full <- get_BIC(logL_full, individuals = nrow(x), clusters = clusters, loci = ncol(x))
        BIC_marginal <- get_BIC(logL_marginal, individuals = nrow(x), clusters = clusters, loci = ncol(x))

        cat(as.character(Sys.time()), ": logL full     = ", logL_full, "\n", sep = "")
        cat(as.character(Sys.time()), ": logL marginal = ", logL_marginal, "\n", sep = "")
        cat(as.character(Sys.time()), ": AIC full      = ", AIC_full, "\n", sep = "")
        cat(as.character(Sys.time()), ": AIC marginal  = ", AIC_marginal, "\n", sep = "")
        cat(as.character(Sys.time()), ": AICc full     = ", AICc_full, "\n", sep = "")
        cat(as.character(Sys.time()), ": AICc marginal = ", AICc_marginal, "\n", sep = "")
        cat(as.character(Sys.time()), ": BIC full      = ", BIC_full, "\n", sep = "")
        cat(as.character(Sys.time()), ": BIC marginal  = ", BIC_marginal, "\n", sep = "")

        logL_full_iterations <- c(logL_full_iterations, logL_full)
        logL_marginal_iterations <- c(logL_marginal_iterations, logL_marginal)
        AIC_full_iterations <- c(AIC_full_iterations, AIC_full)
        AIC_marginal_iterations <- c(AIC_marginal_iterations, AIC_marginal)
        AICc_full_iterations <- c(AICc_full_iterations, AICc_full)
        AICc_marginal_iterations <- c(AICc_marginal_iterations, AICc_marginal)
        BIC_full_iterations <- c(BIC_full_iterations, BIC_full)
        BIC_marginal_iterations <- c(BIC_marginal_iterations, BIC_marginal)
      }

      if (verbose >= 1L) {
        cat(as.character(Sys.time()), ": Iteration ", iter, ", ", sep = "")
        cat("max|vic - vic_old| / max(vic_old) = ", 
          v_max_diff, " / ", last_v_max, " = ", v_gain, " (eps = ", eps, ")\n", sep = "")
      }
      
      if (!is.na(v_gain) && v_gain < eps) {
        converged <- TRUE
        
        if (verbose >= 1L) {
          cat("\n")
          cat(as.character(Sys.time()), ": Stopping after ", iterations_total, 
            " iterations due to convergence, ", eps, " > ", 
            v_gain, "\n\n", sep = "")
        }
        
        break
      }
      
      last_v_max <- max(abs(new_vic_matrix))

      if (converged) {
        break
      }
    }

    ##############################################################################
    # Move centers
    ##############################################################################
    if (verbose >= 1L) {
      cat(as.character(Sys.time()), ": Checking if the central haplotypes, y, are optimal.\n", sep = "")
    }
    new_y <- move_centers(x, y, vic_matrix)
    dist_new_y <- sum(abs(y - new_y))
    
    if (dist_new_y == 0) {
      if (verbose >= 1L) {
        cat(as.character(Sys.time()), ": Central haplotypes, y, were optimal, no need to more EM iterations.\n", sep = "")
      }
      break
    } else {
      if (verbose >= 2L) {      
        cat(as.character(Sys.time()), ": Current central haplotypes, y, not optimal:\n", sep = "")
        print(y)
        cat("New centers:\n")
        print(new_y)
        cat("Differences:\n")
        print(new_y - y)
        cat("Number of stepwise mutations between center configurations = ", dist_new_y, "\n", sep = "")
        cat("Doing EM again with the new centers...\n")
      } else if (verbose >= 1) {
        cat(as.character(Sys.time()), ": Current centers not optimal, moving and making another EM iteration.\n", sep = "")
      }
      
      converged <- FALSE
      changed_center <- c(changed_center, iterations_total)
      centers_iterations[[length(centers_iterations)+1]] <- y
      y <- new_y
      response_vector <- create_response_vector(x, y)
    }
  }

  ################################################################################
  rownames(tau_iterations) <- NULL
  colnames(disclap_parameters) <- colnames(x)
  colnames(y) <- colnames(x)
  ################################################################################
  
  #print(fit$coefficients)
  #print(disclap_parameters)
  #print(str(disclap_parameters))

  if (verbose < 2L) { # If 2L, we already calculated these
    logL_full <- get_loglikelihood_full(fit, clusters, ncol(x), response_vector, weight_vector, tau_norm = tau_vector^(1/ncol(x)))
    logL_marginal <- get_loglikelihood_marginal(x, y, disclap_parameters, tau_vector)
    AIC_full <- get_AIC(logL_full, individuals = nrow(x), clusters = clusters, loci = ncol(x))
    AIC_marginal <- get_AIC(logL_marginal, individuals = nrow(x), clusters = clusters, loci = ncol(x))
    AICc_full <- get_AICc(logL_full, individuals = nrow(x), clusters = clusters, loci = ncol(x))
    AICc_marginal <- get_AICc(logL_marginal, individuals = nrow(x), clusters = clusters, loci = ncol(x))
    BIC_full <- get_BIC(logL_full, individuals = nrow(x), clusters = clusters, loci = ncol(x))
    BIC_marginal <- get_BIC(logL_marginal, individuals = nrow(x), clusters = clusters, loci = ncol(x))
  }
  
  if (clusters > 1L && !converged) {
    msg <- paste("EM did not converge according to the specified eps = ", 
      eps, " (only reached ", v_gain, ")", sep = "")
    warning(msg)
  }

  ans <- list(
    glm_method = glm_method,
#    glm_control_maxit = glm_control_maxit,
#    glm_control_eps = glm_control_eps,
    
    init_y = init_y,
    init_y_method = init_y_method,
    
    #fit = fit,

    converged = converged,
        
    y = y,
    tau = tau_vector,
    v_matrix = vic_matrix,
    disclap_parameters = disclap_parameters,
    glm_coef = fit$coefficients,
    
    model_observations = prod(dim(x)),
    model_parameters = ((clusters * ncol(x)) + (ncol(x) + clusters - 1) + (clusters - 1)),
    iterations = iterations_total,  
    
    logL_full = logL_full,
    logL_marginal = logL_marginal,  
    AIC_full = AIC_full,
    AIC_marginal = AIC_marginal,   
    AICc_full = AICc_full,
    AICc_marginal = AICc_marginal,  
    BIC_full = BIC_full,
    BIC_marginal = BIC_marginal, 

    v_gain_iterations = v_gain_iterations,
    tau_iterations = tau_iterations,
    
#    changed_center = changed_center,
    centers_iterations = centers_iterations,
    
    logL_full_iterations = logL_full_iterations,
    logL_marginal_iterations = logL_marginal_iterations,
    AIC_full_iterations = AIC_full_iterations,
    AIC_marginal_iterations = AIC_marginal_iterations,
    AICc_full_iterations = AICc_full_iterations,
    AICc_marginal_iterations = AICc_marginal_iterations,
    BIC_full_iterations = BIC_full_iterations,
    BIC_marginal_iterations = BIC_marginal_iterations    
  )

  class(ans) <- "disclapmixfit"

  if (verbose >= 1L) {
    cat(as.character(Sys.time()), ": Done.\n", sep = "")
  }  
  
  return(ans)
}

