context("Testing theta2tau and vice versa \n")

theta <- list(beta = c(2, 3, 4), gamma = 2)
theta <- complete_theta(theta)

test_that("theta2tau converts to correct location and scale or mean/variance", {
  tau.tmp <- theta2tau(theta, distname = "t", use.mean.variance = TRUE)
  expect_equivalent(tau.tmp["mu_x"], 2)
  expect_equivalent(tau.tmp["sigma_x"], 3 * sqrt(4 / 2))
 
  tau.tmp.lc <- theta2tau(theta, distname = "t", use.mean.variance = FALSE)
  expect_equivalent(tau.tmp.lc["mu_x"], 2)
  expect_equivalent(tau.tmp.lc["sigma_x"], 3)
})


test_that("tau2theta is inverse of theta2tau", {
  tau.tmp <- theta2tau(theta, distname = "t", use.mean.variance = TRUE)
  
  theta.tmp <- tau2theta(tau.tmp, beta = c(2, 3, 4))
  for (nn in names(theta.tmp)) {
    expect_equivalent(theta.tmp[[nn]], theta[[nn]])
  }
})

test_that("theta2tau computes correct mean/variance for uniform", {
  for (dd in c(0, 0.1, 1)) {
    theta.tmp <- list(beta = c(min = -1, max = 1), delta = dd)
    tau.tmp <- theta2tau(theta.tmp, "unif")
  
    expect_equivalent(tau.tmp["mu_x"], 0)
    expect_equivalent(tau.tmp["sigma_x"], sqrt(1 / 12) * 2)
    expect_equivalent(tau.tmp["delta"], dd)
  }
})