context("Check dgmm.loglik")


# The for proper return of size
n <- 1000
m <- 3
d <- 3
z <- SimulateGMCMData(n = n, m = m, d = d)$z
theta <- rtheta(m = m, d = d)

ans1 <- GMCM:::dgmm.loglik(theta = theta, z = z)
ans2 <- GMCM:::dgmm.loglik(theta = theta, z = z,
                           marginal.loglik = TRUE)

test_that(paste0("dgmm.loglik returns proper size (n = ",
                 n, ", m = ", m, ")"), {
  expect_that(ans1, is_a("matrix"))
  expect_that(nrow(ans1), equals(1))
  expect_that(ncol(ans1), equals(1))
  expect_that(ans2, is_a("matrix"))
  expect_that(nrow(ans2), equals(n))
  expect_that(ncol(ans2), equals(1))
})



# Test wrong inputs
theta$d <- 1e2
expect_that(suppressWarnings(GMCM:::dgmm.loglik(theta = theta, z = z)),
            throws_error("formatted"))

theta$d <- d
expect_that(GMCM:::dgmm.loglik(theta = theta, z = z[,1]),
            throws_error("not a matrix"))
expect_that(GMCM:::dgmm.loglik(theta = theta, z = z[,1:2]),
            throws_error("colums of z does not equal"))
