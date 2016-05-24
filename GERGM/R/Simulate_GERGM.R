# Simulate a gergm
Simulate_GERGM <- function(GERGM_Object,
                           coef = GERGM_Object@theta.par,
                           seed1,
						               possible.stats,
						               verbose = TRUE) {
  # object: an object of class "gergm"

  sample_every <- floor(1/GERGM_Object@thin)
  thetas <- GERGM_Object@theta.par
  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))
  pairs <- t(combn(1:num.nodes, 2))

  # if we are dealing with an undirected network
  undirect_network <- 0
  if (!GERGM_Object@directed_network) {
    undirect_network <- 1
  }

  # if we are dealing with a correlation network
  is_correlation_network <- 0
  if (GERGM_Object@is_correlation_network) {
    is_correlation_network <- 1
    undirect_network <- 1
  }

  # Gibbs Simulation
  if (GERGM_Object@estimation_method == "Gibbs") {
    nets <- Gibbs_Sampler(GERGM_Object,
                          thetas,
                          MCMC.burnin = GERGM_Object@burnin,
                          num.draws = GERGM_Object@number_of_simulations,
                          thin = GERGM_Object@thin,
                          start = NULL,
                          num.nodes = num.nodes,
                          directed = TRUE,
                          possible.stats = possible.stats)
    # Calculate the network statistics over all of the simulated networks
    h.statistics <- t(apply(nets, 3, h2,
                            triples = triples,
                            statistics = rep(1, length(possible.stats)),
                            alphas = rep(1, length(possible.stats)),
                            together = GERGM_Object@downweight_statistics_together))
    acceptance.rate <- NULL
  }

  # Metropolis Hastings Simulation
  if (GERGM_Object@estimation_method == "Metropolis") {
    #need to put the thetas into a full length vector for MH function
    stat.indx <- which(GERGM_Object@stats_to_use > 0)
    #cat("stat.idx",stat.indx,"\n" )
    full_thetas <- rep(0, length(GERGM_Object@stats_to_use))
    for (i in 1:length(thetas)) {
      full_thetas[stat.indx[i]] <- thetas[i]
    }

    #cat("Current Theta Estimates:",thetas,"\n")
    store <- ceiling((GERGM_Object@number_of_simulations + GERGM_Object@burnin)/sample_every)
    nsim <- GERGM_Object@number_of_simulations + GERGM_Object@burnin
    dw <- as.numeric(GERGM_Object@downweight_statistics_together)
    samples <- Metropolis_Hastings_Sampler(
      number_of_iterations = nsim,
      shape_parameter = GERGM_Object@proposal_variance,
      number_of_nodes = num.nodes,
      statistics_to_use = GERGM_Object@stats_to_use,
      initial_network = GERGM_Object@bounded.network,
      take_sample_every = sample_every,
      thetas = full_thetas,
      triples = triples - 1,
      pairs = pairs - 1,
      alphas = GERGM_Object@weights,
      together = dw,
      seed = seed1,
      number_of_samples_to_store = store,
      using_correlation_network = is_correlation_network,
      undirect_network = undirect_network)
    # keep only the networks after the burnin
    start <- floor(GERGM_Object@burnin/sample_every) + 1
    end <- length(samples[[3]][,1])
    nets <- samples[[2]][, , start:end]
    # Note: these statistics will be the adjusted statistics (for use in the
    # MCMCMLE procedure)

    # more markov chain diagnostics
    average_log_prob_accept <- mean(samples[[5]])
    cat("Average log probability of accepting a proposal:",
        average_log_prob_accept,
        ".\nStandard deviation of log probability of accepting proposal:",
        sd(samples[[5]]),"\n")
    if (!is.finite(average_log_prob_accept)) {
      warning("It appears there is a problem with Metropolis Hastings, consider increasing proposal variance.")
    }

    GERGM_Object <- store_console_output(GERGM_Object,
      paste("Average log probability of accepting a proposal:",
            average_log_prob_accept,
            ".\nStandard deviation of log probability of accepting proposal:",
            sd(samples[[5]]),"\n"))

    average_edge_weight <- mean(samples[[4]])
    cat("Average (constrained) simulated network density:",
        average_edge_weight, "\n")
    GERGM_Object <- store_console_output(GERGM_Object,
      paste("Average (constrained) simulated network density:",
            average_edge_weight, "\n"))


    h.statistics <- samples[[3]][start:end,]
    acceptance.rate <- mean(samples[[1]])
    if (verbose) {
      cat("Metropolis Hastings Acceptance Rate (target = ",
          GERGM_Object@target_accept_rate," ): ",
          acceptance.rate, "\n", sep = "")
    }
    GERGM_Object <- store_console_output(GERGM_Object,
      paste("Metropolis Hastings Acceptance Rate (target = ",
            GERGM_Object@target_accept_rate,"):",
            acceptance.rate, "\n", sep = ""))

  }
  h.statistics = data.frame(out2stars = h.statistics[, 1],
                            in2stars = h.statistics[, 2],
                            ctriads = h.statistics[, 3],
                            mutual = h.statistics[, 4],
                            ttriads = h.statistics[, 5],
                            edges = h.statistics[, 6])

  GERGM_Object@MCMC_output = list(Networks = nets,
                            Statistics = h.statistics,
                            Acceptance.rate = acceptance.rate)
  return(GERGM_Object)
}
