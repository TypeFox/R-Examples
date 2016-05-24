hysteresis_parallel <- function(i,
                                GERGM_Object,
                                initial_density,
                                possible_structural_terms,
                                seed,
                                steps,
                                observed_density,
                                range){

  # figure out the range of values for each parameter
  current_theta <- GERGM_Object@theta.par[i]
  theta_se <- GERGM_Object@theta.coef[2,i]
  min_val <- current_theta - range * theta_se
  max_val <- current_theta + range * theta_se
  hysteresis_values <- seq(min_val, max_val, length.out = 2 * steps + 1)
  last_network <- floor(GERGM_Object@number_of_simulations*GERGM_Object@thin)
  network_densities <- matrix(0,
    nrow = ceiling(GERGM_Object@number_of_simulations*GERGM_Object@thin),
    ncol = 2*length(hysteresis_values))
  #GERGM_Object@bounded.network <- bounded_network
  n_nodes <- nrow(GERGM_Object@bounded.network)
  zero_net <- matrix(initial_density,n_nodes,n_nodes)
  GERGM_Object@bounded.network <- zero_net
  which_term <- which(GERGM_Object@stats_to_use > 0)[i]
  cur_term <- possible_structural_terms[which_term]
  # tell the user what is going on
  cat("Currently simulating networks while varying the",
      cur_term,"parameter from:",min_val,"to",max_val,"for a total of",
      length(hysteresis_values),"simulations...\n")

  # loop over values for theta
  column_counter <- 1
  for(j in 1:length(hysteresis_values)){

    # set the current value
    GERGM_Object@theta.par[i] <- hysteresis_values[j]
    cat("Current theta values:",GERGM_Object@theta.par,"\n")
    # simulate networks
    GERGM_Object <- Simulate_GERGM(
      GERGM_Object,
      seed1 = seed,
      possible.stats = possible_structural_terms)

    # assign the last simulated network as the starting network for the next
    # simulation
    GERGM_Object@bounded.network <- GERGM_Object@MCMC_output$Networks[,,last_network]
    # save the densities
    nr <- nrow(GERGM_Object@network)
    normalizer <- nr * (nr - 1)
    network_densities[,column_counter] <- GERGM_Object@MCMC_output$Statistics$edges/normalizer
    column_counter <- column_counter + 1
  }
  # now back down
  for (j in length(hysteresis_values):1) {

    # set the current value
    GERGM_Object@theta.par[i] <- hysteresis_values[j]
    cat("Current theta values:",GERGM_Object@theta.par,"\n")
    # simulate networks
    GERGM_Object <- Simulate_GERGM(
      GERGM_Object,
      seed1 = seed,
      possible.stats = possible_structural_terms)

    # assign the last simulated network as the starting network for the next
    # simulation
    GERGM_Object@bounded.network <- GERGM_Object@MCMC_output$Networks[,,last_network]

    # save the densities
    nr <- nrow(GERGM_Object@network)
    normalizer <- nr * (nr - 1)
    network_densities[,column_counter] <- GERGM_Object@MCMC_output$Statistics$edges/normalizer
    column_counter <- column_counter + 1
  }

  # reset the theta value
  GERGM_Object@theta.par[i] <- current_theta
  mean_densities <- apply(network_densities,2,mean)
  thetas <- c(hysteresis_values, rev(hysteresis_values))

  hysteresis_dataframe <- data.frame(theta_values = thetas,
                                     mean_densities = mean_densities)

  Hysteresis_Results <- list(network_densities = network_densities,
                                  mean_densities = mean_densities,
                                  theta_values = hysteresis_values,
                                  hysteresis_dataframe = hysteresis_dataframe,
                                  observed_density = observed_density,
                                  term = cur_term)

  return(Hysteresis_Results)
}


#intitalizes snowfall session
#     snowfall::sfInit(parallel = TRUE, cpus = cores)
#
#     #check to see if we are running in parallel
#     if(snowfall::sfParallel())
#       cat( "Running in parallel mode on", snowfall::sfCpus(), "nodes.\n" )
#     else
#       cat( "Running in sequential mode.\n" )
#
#     #export all packages and libraries currently loaded in workspace
#     for (i in 1:length(.packages())){
#       eval(call("sfLibrary", (.packages()[i]), character.only = TRUE))
#     }
#
#     # apply our problem across the cluster using hte indexes we have determined
#     # and load balancing
#     # Export a list of R data objects
#     snowfall::sfExport("GERGM_Object",
#                        "initial_density",
#                        "possible_structural_terms",
#                        "networks_to_simulate",
#                        "simulation_method",
#                        "burnin",
#                        "thin",
#                        "proposal_variance",
#                        "seed",
#                        "steps",
#                        "observed_density",
#                        "range")
#                        print(vec)
#     results <- snowfall::sfClusterApplyLB(x = vec,
#       fun = hysteresis_parallel)

#stop the cluster when we are done -- this is very important and must be
#done manually every time
# snowfall::sfStop()

