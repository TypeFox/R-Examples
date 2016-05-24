test_that("Hyperparameter optimization works", {
  skip_on_cran()

  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)

  formula <- net ~ mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodematch("Type",base = "A")

  test <- gergm(formula,
                covariate_data = node_level_covariates,
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

  check_against <- c(1.288, -0.075, -0.017, -0.025,  3.127, 0.132, -1.837)
  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  expect_equal(check, check_against)

})





test_that("Model works for correlation networks", {
  skip_on_cran()
  skip("Skipping test as it can only be run in the global environment.")

  set.seed(12345)
  #Function to generating a random positive-definite matrix with user-specified positive
  #eigenvalues
  # If eigenvalues are not specified, they are generated from a uniform
  #distribution

  Posdef <- function (n, ev = runif(n, 0, 10))
  {
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
  }

  #Generate a random correlation matrix of dimension 10 x 10
  x <- rnorm(10)

  pdmat <- Posdef(n = 10)
  correlations <- pdmat / max(abs(pdmat))
  diag(correlations) <- 1
  net <- (correlations + t(correlations)) / 2
  colnames(net) <- rownames(net) <- letters[1:10]

  formula <- net ~ edges + ttriads

  test <- gergm(formula,
                normalization_type = "division",
                network_is_directed = FALSE,
                use_MPLE_only = FALSE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 40000,
                thin = 1/10,
                proposal_variance = 0.03,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 10000,
                seed = 456,
                convergence_tolerance = 0.01,
                MPLE_gain_factor = 0,
                force_x_theta_updates = 4,
                force_x_lambda_updates = 2,
                using_correlation_network = TRUE,
                omit_intercept_term = TRUE)

})
