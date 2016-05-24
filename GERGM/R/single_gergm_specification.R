# helper function to run GERGM's in parallel
single_gergm_specification <- function(i,
                                       num_specifications,
                                       formula_list,
                                       observed_network_list,
                                       covariate_data_list,
                                       network_data_list,
                                       normalization_type,
                                       network_is_directed,
                                       use_MPLE_only,
                                       transformation_type,
                                       estimation_method,
                                       maximum_number_of_lambda_updates,
                                       maximum_number_of_theta_updates,
                                       number_of_networks_to_simulate,
                                       thin,
                                       proposal_variance,
                                       downweight_statistics_together,
                                       MCMC_burnin,
                                       seed,
                                       convergence_tolerance,
                                       MPLE_gain_factor,
                                       acceptable_fit_p_value_threshold,
                                       force_x_theta_updates,
                                       force_x_lambda_updates,
                                       output_directory,
                                       output_name,
                                       generate_plots,
                                       verbose,
                                       omit_intercept_term,
                                       hyperparameter_optimization ,
                                       target_accept_rate,
                                       ...){

  # 0. go through and assign all variables before calling GERGM

  # get the formula
  if (length(formula_list) == num_specifications) {
    formula <- as.formula(formula_list[[i]])
  } else if (length(formula_list) == 1) {
    if (class(formula_list) == "list") {
      formula <- as.formula(formula_list[[1]])
    } else if (class(formula_list) == "character" |
               class(formula_list) == "formula") {
      formula <- as.formula(formula_list)
    } else {
      stop("formula_list must be provided as either a list of formulas, a string,
           or a formula")
    }
  } else {
    stop(paste("formula_list must either be of length 1 or of length",
               num_specifications))
  }


  # figure out what the user called the dependent network variable and assign the
  # matrix in this list to that value so things work.
  if (class(formula) != "formula") {
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])

  if(is.null(observed_network_list)){
    net <- dynGet(as.character(lhs),
                  ifnotfound = get(as.character(lhs)))
    assign(lhs,net)
  }else{
    if (class(observed_network_list) == "matrix") {
      assign(lhs,observed_network_list)
    } else {
      if (length(observed_network_list) == num_specifications) {
        assign(lhs,observed_network_list[[i]])
      } else if (length(observed_network_list) == 1) {
        if (class(observed_network_list) == "list") {
          assign(lhs,observed_network_list[[1]])
        } else {
          stop("observed_network_list must be provided as either a list of
             matrices or a single numeric matrix")
        }
      } else {
        stop(paste("observed_network_list must either be of length 1 or of length",
                   num_specifications))
      }
    }
  }

  # get the covariate data frame
  if (is.null(covariate_data_list)) {
    covariate_data <- NULL
  } else {
    if (class(covariate_data_list) == "data.frame") {
      covariate_data <- covariate_data_list
    } else {
      if (length(covariate_data_list) == num_specifications) {
        covariate_data <- covariate_data_list[[i]]
      } else if (length(covariate_data_list) == 1) {
        if (class(covariate_data_list) == "list") {
          covariate_data <- covariate_data_list[[1]]
        } else {
          stop("covariate_data_list must be provided as either a list of
               data.frames or a single data.frame")
        }
        } else {
          stop(paste("covariate_data_list must either be of length 1 or of length",
                     num_specifications))
      }
    }
  }

  if (!is.null(network_data_list)) {
    # assign all network covariate data objects
    if (length(network_data_list) == num_specifications) {
      temp <- network_data_list[[i]]
      for(j in 1:length(temp)){
        assign(names(temp)[j],temp[[j]])
      }
    } else if (length(network_data_list) == 1) {
      if (class(network_data_list) == "list" &
          class(network_data_list[[1]]) == "matrix") {
        for(j in 1:length(network_data_list)){
          assign(names(network_data_list)[j],network_data_list[[j]])
        }
      } else if (class(network_data_list) == "list") {
        temp <- network_data_list[[1]]
        for(j in 1:length(temp)){
          assign(names(temp)[j],temp[[j]])
        }
      } else if (class(network_data_list) == "symbol") {
        net <- dynGet(as.character(network_data_list),
                      ifnotfound = get(as.character(network_data_list)))
        assign(as.character(network_data_list),net)
      } else {
        stop("network_data_list must be provided as either a list of
             lists or a single matrix object")
      }
    } else {
      for(j in 1:length(temp)){
        assign(names(network_data_list)[j],network_data_list[[j]])
      }
    }
  }

  # 1. now we deal with all of the variables that can be either fixed or vary
  # across specifications. Note that we leave on the case where length() = 1,
  #  implictly passing single arguments through.

  if (length(normalization_type) == num_specifications) {
    if (class(normalization_type) == "list") {
      normalization_type <- normalization_type[[i]]
    } else if (class(normalization_type) == "character") {
      normalization_type <- normalization_type[i]
    } else {
      stop("normalization_type must either be a character vector or list of strings.")
    }
  } else if (length(normalization_type) > 2) {
    stop("normalization_type must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(network_is_directed) == num_specifications) {
    if (class(network_is_directed) == "list") {
      normalization_type <- network_is_directed[[i]]
    } else if (class(network_is_directed) == "logical") {
      network_is_directed <- network_is_directed[i]
    } else {
      stop("network_is_directed must either be a logical vector or list of logicals.")
    }
  } else if (length(network_is_directed) != 1) {
    stop("network_is_directed must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(use_MPLE_only) == num_specifications) {
    if (class(use_MPLE_only) == "list") {
      normalization_type <- use_MPLE_only[[i]]
    } else if (class(use_MPLE_only) == "logical") {
      use_MPLE_only <- use_MPLE_only[i]
    } else {
      stop("use_MPLE_only must either be a logical vector or list of logicals.")
    }
  } else if (length(use_MPLE_only) != 1) {
    stop("use_MPLE_only must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(transformation_type) == num_specifications) {
    if (class(transformation_type) == "list") {
      normalization_type <- transformation_type[[i]]
    } else if (class(transformation_type) == "character") {
      transformation_type <- transformation_type[i]
    } else {
      stop("transformation_type must either be a character vector or list of strings.")
    }
  } else if (length(transformation_type) > 4) {
    stop("transformation_type must either be the same length as the number of specifications or of length four, in which case it will be the same across all specifications.")
  }

  if (length(estimation_method) == num_specifications) {
    if (class(estimation_method) == "list") {
      normalization_type <- estimation_method[[i]]
    } else if (class(estimation_method) == "character") {
      estimation_method <- estimation_method[i]
    } else {
      stop("estimation_method must either be a character vector or list of strings.")
    }
  } else if (length(estimation_method) > 2) {
    # pass through default
    stop("estimation_method must either be the same length as the number of specifications or of length two, in which case it will be the same across all specifications.")
  }

  if (length(maximum_number_of_lambda_updates) == num_specifications) {
    if (class(maximum_number_of_lambda_updates) == "list") {
      normalization_type <- maximum_number_of_lambda_updates[[i]]
    } else if (class(maximum_number_of_lambda_updates) == "numeric") {
      maximum_number_of_lambda_updates <- maximum_number_of_lambda_updates[i]
    } else {
      stop("maximum_number_of_lambda_updates must either be a numeric vector or list of numbers.")
    }
  } else if (length(maximum_number_of_lambda_updates) != 1) {
    stop("maximum_number_of_lambda_updates must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(maximum_number_of_theta_updates) == num_specifications) {
    if (class(maximum_number_of_theta_updates) == "list") {
      normalization_type <- maximum_number_of_theta_updates[[i]]
    } else if (class(maximum_number_of_theta_updates) == "numeric") {
      maximum_number_of_theta_updates <- maximum_number_of_theta_updates[i]
    } else {
      stop("maximum_number_of_theta_updates must either be a numeric vector or list of numbers.")
    }
  } else if (length(maximum_number_of_theta_updates) != 1) {
    stop("maximum_number_of_theta_updates must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(number_of_networks_to_simulate) == num_specifications) {
    if (class(number_of_networks_to_simulate) == "list") {
      normalization_type <- number_of_networks_to_simulate[[i]]
    } else if (class(number_of_networks_to_simulate) == "numeric") {
      number_of_networks_to_simulate <- number_of_networks_to_simulate[i]
    } else {
      stop("number_of_networks_to_simulate must either be a numeric vector or list of numbers.")
    }
  } else if (length(number_of_networks_to_simulate) != 1) {
    stop("number_of_networks_to_simulate must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(thin) == num_specifications) {
    if (class(thin) == "list") {
      normalization_type <- thin[[i]]
    } else if (class(thin) == "numeric") {
      thin <- thin[i]
    } else {
      stop("thin must either be a numeric vector or list of numbers.")
    }
  } else if (length(thin) != 1) {
    stop("thin must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(proposal_variance) == num_specifications) {
    if (class(proposal_variance) == "list") {
      normalization_type <- proposal_variance[[i]]
    } else if (class(proposal_variance) == "numeric") {
      proposal_variance <- proposal_variance[i]
    } else {
      stop("proposal_variance must either be a numeric vector or list of numbers.")
    }
  } else if (length(proposal_variance) != 1) {
    stop("proposal_variance must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(downweight_statistics_together) == num_specifications) {
    if (class(downweight_statistics_together) == "list") {
      normalization_type <- downweight_statistics_together[[i]]
    } else if (class(downweight_statistics_together) == "logical") {
      downweight_statistics_together <- downweight_statistics_together[i]
    } else {
      stop("downweight_statistics_together must either be a logical vector or list of logicals.")
    }
  } else if (length(downweight_statistics_together) != 1) {
    stop("downweight_statistics_together must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(MCMC_burnin) == num_specifications) {
    if (class(MCMC_burnin) == "list") {
      normalization_type <- MCMC_burnin[[i]]
    } else if (class(MCMC_burnin) == "numeric") {
      MCMC_burnin <- MCMC_burnin[i]
    } else {
      stop("MCMC_burnin must either be a numeric vector or list of numbers.")
    }
  } else if (length(MCMC_burnin) != 1) {
    stop("MCMC_burnin must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(seed) == num_specifications) {
    if (class(seed) == "list") {
      normalization_type <- seed[[i]]
    } else if (class(seed) == "numeric") {
      seed <- seed[i]
    } else {
      stop("seed must either be a numeric vector or list of numbers.")
    }
  } else if (length(seed) != 1) {
    stop("seed must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(convergence_tolerance) == num_specifications) {
    if (class(convergence_tolerance) == "list") {
      normalization_type <- convergence_tolerance[[i]]
    } else if (class(convergence_tolerance) == "numeric") {
      convergence_tolerance <- convergence_tolerance[i]
    } else {
      stop("convergence_tolerance must either be a numeric vector or list of numbers.")
    }
  } else if (length(convergence_tolerance) != 1) {
    stop("convergence_tolerance must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(MPLE_gain_factor) == num_specifications) {
    if (class(MPLE_gain_factor) == "list") {
      normalization_type <- MPLE_gain_factor[[i]]
    } else if (class(MPLE_gain_factor) == "numeric") {
      MPLE_gain_factor <- MPLE_gain_factor[i]
    } else {
      stop("MPLE_gain_factor must either be a numeric vector or list of numbers.")
    }
  } else if (length(MPLE_gain_factor) != 1) {
    stop("MPLE_gain_factor must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(acceptable_fit_p_value_threshold) == num_specifications) {
    if (class(acceptable_fit_p_value_threshold) == "list") {
      normalization_type <- acceptable_fit_p_value_threshold[[i]]
    } else if (class(acceptable_fit_p_value_threshold) == "numeric") {
      acceptable_fit_p_value_threshold <- acceptable_fit_p_value_threshold[i]
    } else {
      stop("acceptable_fit_p_value_threshold must either be a numeric vector or list of numbers.")
    }
  } else if (length(acceptable_fit_p_value_threshold) != 1) {
    stop("acceptable_fit_p_value_threshold must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(force_x_theta_updates) == num_specifications) {
    if (class(force_x_theta_updates) == "list") {
      normalization_type <- force_x_theta_updates[[i]]
    } else if (class(force_x_theta_updates) == "numeric") {
      force_x_theta_updates <- force_x_theta_updates[i]
    } else {
      stop("force_x_theta_updates must either be a numeric vector or list of numbers.")
    }
  } else if (length(force_x_theta_updates) != 1) {
    stop("force_x_theta_updates must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(force_x_lambda_updates) == num_specifications) {
    if (class(force_x_lambda_updates) == "list") {
      normalization_type <- force_x_lambda_updates[[i]]
    } else if (class(force_x_lambda_updates) == "numeric") {
      force_x_lambda_updates <- force_x_lambda_updates[i]
    } else {
      stop("force_x_lambda_updates must either be a numeric vector or list of numbers.")
    }
  } else if (length(force_x_lambda_updates) != 1) {
    stop("force_x_lambda_updates must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(verbose) == num_specifications) {
    if (class(verbose) == "list") {
      normalization_type <- verbose[[i]]
    } else if (class(verbose) == "logical") {
      verbose <- verbose[i]
    } else {
      stop("verbose must either be a logical vector or list of logicals.")
    }
  } else if (length(verbose) != 1) {
    stop("verbose must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(omit_intercept_term) == num_specifications) {
    if (class(omit_intercept_term) == "list") {
      normalization_type <- omit_intercept_term[[i]]
    } else if (class(omit_intercept_term) == "logical") {
      omit_intercept_term <- omit_intercept_term[i]
    } else {
      stop("omit_intercept_term must either be a logical vector or list of logicals.")
    }
  } else if (length(omit_intercept_term) != 1) {
    stop("omit_intercept_term must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(hyperparameter_optimization) == num_specifications) {
    if (class(hyperparameter_optimization) == "list") {
      normalization_type <- hyperparameter_optimization[[i]]
    } else if (class(hyperparameter_optimization) == "logical") {
      hyperparameter_optimization <- hyperparameter_optimization[i]
    } else {
      stop("hyperparameter_optimization must either be a logical vector or list of logicals.")
    }
  } else if (length(hyperparameter_optimization) != 1) {
    stop("hyperparameter_optimization must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  if (length(target_accept_rate) == num_specifications) {
    if (class(target_accept_rate) == "list") {
      normalization_type <- target_accept_rate[[i]]
    } else if (class(target_accept_rate) == "numeric") {
      target_accept_rate <- target_accept_rate[i]
    } else {
      stop("target_accept_rate must either be a numeric vector or list of numbers.")
    }
  } else if (length(target_accept_rate) != 1) {
    stop("target_accept_rate must either be the same length as the number of specifications or of length one, in which case it will be the same across all specifications.")
  }

  Result <- gergm(formula = formula,
    covariate_data = covariate_data,
    normalization_type = normalization_type,
    network_is_directed = network_is_directed,
    use_MPLE_only = use_MPLE_only,
    transformation_type = transformation_type,
    estimation_method = estimation_method,
    maximum_number_of_lambda_updates = maximum_number_of_lambda_updates,
    maximum_number_of_theta_updates = maximum_number_of_theta_updates,
    number_of_networks_to_simulate = number_of_networks_to_simulate,
    thin = thin,
    proposal_variance = proposal_variance,
    downweight_statistics_together = downweight_statistics_together,
    MCMC_burnin = MCMC_burnin,
    seed = seed,
    convergence_tolerance = convergence_tolerance,
    MPLE_gain_factor = MPLE_gain_factor,
    acceptable_fit_p_value_threshold = acceptable_fit_p_value_threshold,
    force_x_theta_updates = force_x_theta_updates,
    force_x_lambda_updates = force_x_lambda_updates,
    output_directory = output_directory,
    output_name = output_name,
    generate_plots = generate_plots,
    verbose = verbose,
    omit_intercept_term = omit_intercept_term,
    hyperparameter_optimization = hyperparameter_optimization,
    target_accept_rate = target_accept_rate,
    ... = ...)

  return(Result)
}
