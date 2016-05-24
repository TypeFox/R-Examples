context("Check meta.gmcm.loglik")

n <- 1000
tpar <- c(tpie1 = rnorm(1, sd = 2), tmu = rchisq(1, 5),
          tsigma = rchisq(1, 4),    trho = rnorm(1, sd = 4))

u <- SimulateGMCMData(n = n, par = GMCM:::tt(tpar, d = 3, positive.rho = TRUE),
                      d = 3, rescale = TRUE, positive.rho = TRUE)$u

ans <- GMCM:::meta.gmcm.loglik(tpar = tpar, u = u)
ans2 <-
  GMCM:::meta.gmcm.loglik(tpar = GMCM:::tt(tpar, d = 3, positive.rho = TRUE),
                          u = u, rescale = FALSE)

test_that("Check meta.gmcm.loglik returns proper formatted output", {
  expect_that(is.numeric(ans),  is_true())
  expect_that(length(ans), equals(1))
  expect_that(ans, equals(ans2))
})

# Test degenerate input
