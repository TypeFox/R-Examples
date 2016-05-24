context("Testing support of distribution\n")

test_that("support is correct for corner cases", {
  for (gg in c(0, 0.1)) {
    expect_equal(get_support(tau = c(mu_x = 0, sigma_x = 1, gamma = gg), TRUE),
                 c(0, Inf),
                 check.names = FALSE)
  }
  
  for (dd in c(0, 0.2, 1)) {
    expect_equal(get_support(tau = c(mu_x = 0, sigma_x = 1, delta = dd), TRUE),
                 c(0, Inf),
                 check.names = FALSE)
  }
  
  for (dd in c(0, 0.2, 1)) {
    expect_equal(get_support(tau = c(mu_x = 0, sigma_x = 1, delta = dd), FALSE),
                 c(-Inf, Inf),
                 check.names = FALSE)
  }
})

test_that("get_support is correct for uniform distribution", {

  theta.tmp <- list(beta = c(min = -1, max = 1), delta = 0)
  tau.tmp <- theta2tau(theta.tmp, distname = "unif")
  expect_identical(get_support(tau = tau.tmp, input.bounds = c(-1, 1)),
                   c(lower = -1, upper = 1))
})