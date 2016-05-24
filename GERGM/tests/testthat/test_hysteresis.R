test_that("Simple model with no covariates runs", {
  skip_on_cran()
  ########################### 1. No Covariates #############################
  # Preparing an unbounded network without covariates for gergm estimation #
  skip("Skipping test as it can only be run in the global environment.")

  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]

  # three parameter model
  formula <- net ~  mutual +  ttriads

  test <- gergm(formula,
                normalization_type = "division",
                network_is_directed = TRUE,
                use_MPLE_only = FALSE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 40000,
                thin = 1/100,
                proposal_variance = 0.5,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 10000,
                seed = 456,
                convergence_tolerance = 0.01,
                MPLE_gain_factor = 0,
                force_x_theta_updates = 4)

  test2 <- hysteresis(test,
                      networks_to_simulate = 20000,
                      burnin = 10000,
                      range = 2,
                      steps = 30,
                      initial_density = 0.2,
                      simulation_method = "Metropolis",
                      proposal_variance = 0.5,
                      seed = 12345,
                      thin = 1/10)

  test3 <- hysteresis(test,
                      networks_to_simulate = 20000,
                      burnin = 10000,
                      range = 2,
                      steps = 30,
                      initial_density = 0.2,
                      simulation_method = "Metropolis",
                      proposal_variance = 0.5,
                      seed = 12345,
                      thin = 1/10,
                      parallel = T)

  hysteresis_plot(test3)

})



