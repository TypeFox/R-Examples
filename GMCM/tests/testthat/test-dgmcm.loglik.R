context("Check dgmcm.loglik")


# The for proper return of size
for (d in 2:4) {
  for (n in c(100, 1232, 2313)) {
    for (m in 2:5) {

      suppressWarnings({
        u <- SimulateGMCMData(n = n, m = m, d = d)$u
      })
      theta <- rtheta(m = m, d = d)


      ans1 <- GMCM:::dgmcm.loglik(theta = theta, u = u)
      ans2 <- GMCM:::dgmcm.loglik(theta = theta, u = u, marginal.loglik = TRUE)

      test_that(paste0("dgmcm.loglik returns proper size (n = ",
                       n, ", m = ", m, ")"), {
        expect_that(ans1, is_a("numeric"))
        expect_that(length(ans1), equals(1L))
        expect_that(ans2, is_a("matrix"))
        expect_that(length(ans2), equals(n))
        expect_that(nrow(ans2), equals(n))
        expect_that(ncol(ans2), equals(1))
      })
    }
  }
}


