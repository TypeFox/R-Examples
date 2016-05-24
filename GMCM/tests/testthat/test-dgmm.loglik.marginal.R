context("Check dgmm.loglik.marginal")


# The for proper return of size
n <- 1000
m <- 3
d <- 3
z <- SimulateGMCMData(n = n, m = m, d = d)$z
theta <- rtheta(m = m, d = d)

ans1 <- GMCM:::dgmm.loglik.marginal(theta = theta, x = z)
ans2 <- GMCM:::dgmm.loglik.marginal(theta = theta, x = z,
                                    marginal.loglik = FALSE)

test_that(paste0("dgmm.loglik.marginal returns proper size (n = ",
                 n, ", m = ", m, ")"), {
  expect_that(ans1, is_a("matrix"))
  expect_that(nrow(ans1), equals(n))
  expect_that(ncol(ans1), equals(d))
  expect_that(ans2, is_a("matrix"))
  expect_that(nrow(ans2), equals(1))
  expect_that(ncol(ans2), equals(d))
})



# Test wrong inputs
theta$d <- 1e2
expect_that(suppressWarnings(GMCM:::dgmm.loglik.marginal(theta = theta, x = z)),
            throws_error("formatted"))

theta$d <- d
expect_that(GMCM:::dgmm.loglik.marginal(theta = theta, x = z[,1]),
            throws_error("not a matrix"))
expect_that(GMCM:::dgmm.loglik.marginal(theta = theta, x = z[,1:2]),
            throws_error("colums of x does not equal"))
