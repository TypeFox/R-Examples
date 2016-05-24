context("Testing Gaussianize \n")

test_that("get_input and get_output are inverse of each other", {
  tau.tmp <- c(mu_x = 0, sigma_x = 1, delta = 0.2)
  xx <- rnorm(100)
  yy <- get_output(xx, tau.tmp)
  xx.hat <- get_input(yy, tau.tmp)
  # they must be equal (except for numerical issues)
  expect_equal(lp_norm(xx - xx.hat, 1), 0) 
})

yy <- rLambertW(n = 100, distname = "normal", theta = list(beta = c(3, 4), gamma = 0.2))
fit.gmm <- IGMM(yy, type = "s")

test_that("return.u actually returns a the zero-mean, unit-variance vector", {
  
  xx.input <- get_input(fit.gmm)
  xx.input.direct <- get_input(yy, tau = fit.gmm$tau)
  expect_equal(xx.input, xx.input.direct)
  
  input.xx.uu <- get_input(fit.gmm, return.u = TRUE)
  expect_true(inherits(input.xx.uu, "list"))
  expect_identical(names(input.xx.uu), c("u", "x"))
  
  expect_equal(cor(input.xx.uu$u, input.xx.uu$x), 1)
  expect_equal(mean(input.xx.uu$u), 0, tol = 1e-4)
  expect_equal(sd(input.xx.uu$u), 1, tol = 1e-4)
})