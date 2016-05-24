context("Check EStep")

set.seed(123021930)
for (n in c(100, 200, 300, 400)) {
  for (m in 1:4) {

    sim <- GMCM:::SimulateGMMData(n = n, m = m)
    init.theta <- GMCM:::rtheta(m = m)  # Generate starting parameters
    es <- GMCM:::EStep(sim$z, init.theta)

    test_that("EStep returns proper formatted output", {
      expect_that(is.matrix(es),  is_true())
      expect_that(is.numeric(es), is_true())
      expect_that(nrow(es), equals(n))
      expect_that(ncol(es), equals(m))
    })
  }
}

# Test more parameters!
