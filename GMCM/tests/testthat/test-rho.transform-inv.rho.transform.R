context("Check rho.transform and inv.rho.transform")

d <- 4
n <- 1000
rho  <- runif(n, -1/(d-1), 1)
a    <- GMCM:::rho.transform(rho, d)
rho2 <- GMCM:::inv.rho.transform(a, d)

test_that("rho.transform and inv.rho.transform returns proper formatted output", {
  expect_that(is.numeric(a),  is_true())
  expect_that(length(a),  equals(n))
  expect_that(is.numeric(rho2),  is_true())
  expect_that(length(rho2),  equals(n))
  expect_that(rho - rho2,  equals(rep(0, n)))
})

# Test degenerate input
