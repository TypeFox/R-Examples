test_that("Parallel GERGMs", {
  skip_on_cran()

  skip("Tests should not call parallel processes")

  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)

  network_data_list <- list(network_covariate = network_covariate)

  formula <- net ~ mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodematch("Type",base = "A")
  formula2 <- net ~ mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodemix("Type",base = "A")

  form_list <- list(f1 = formula,
                    f2 = formula2)

  testl <- parallel_gergm(formula_list = form_list,
                observed_network_list = net,
                covariate_data_list = node_level_covariates,
                network_data_list = network_data_list,
                cores = 2,
                network_is_directed = TRUE,
                use_MPLE_only = FALSE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 100000,
                thin = 1/100,
                proposal_variance = 0.1,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 50000,
                seed = 456,
                convergence_tolerance = 0.01,
                MPLE_gain_factor = 0,
                force_x_theta_updates = 2,
                hyperparameter_optimization = TRUE
  )

  test <- testl[[1]]
  check_against <- c(1.288, -0.075, -0.017, -0.025,  3.127, 0.132, -1.837)
  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  expect_equal(check, check_against)

  test <- testl[[2]]
  check_against <- c(0.835, -0.073, -0.016, -0.026, -0.024, -0.056, -0.055,
                     -0.035, 0.002, -0.040, -0.050,  3.056,  0.129, -1.931)
  expect_equal(c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3)), check_against)

})

