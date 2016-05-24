context("Check SimulateGMCMData and SimualteGMMData")

for (d in 2:4) {
  for (n in c(10, 100, 1000)) {
    suppressWarnings({
      gmm <- SimulateGMMData(n = n, theta = rtheta(d = d))
    })

    test_that("SimulateGMMData returns proper format", {
      expect_that(is.list(gmm), is_true())
      expect_that(length(gmm), equals(3))
      expect_that(names(gmm), equals(c("z", "K", "theta")))
      # Check z
      expect_that(is.matrix(gmm$z), is_true())
      expect_that(is.numeric(gmm$z), is_true())
      expect_that(dim(gmm$z), equals(c(n, d)))
      # Check K
      expect_that(is.numeric(gmm$K), is_true())
      expect_that(length(gmm$K), equals(n))
      # Check theta
      expect_that(is.theta(gmm$theta), is_true())
    })

    suppressWarnings({
      gmcm <- SimulateGMCMData(n = n, theta = rtheta(d = d))
    })

    test_that("SimulateGMCMData returns proper format", {
      expect_that(is.list(gmcm), is_true())
      expect_that(length(gmcm), equals(4))
      expect_that(names(gmcm), equals(c("u", "z", "K", "theta")))
      # Check u
      expect_that(is.matrix(gmcm$u), is_true())
      expect_that(is.numeric(gmcm$u), is_true())
      expect_that(dim(gmcm$u), equals(c(n, d)))
      # Check z
      expect_that(is.matrix(gmcm$z), is_true())
      expect_that(is.numeric(gmcm$z), is_true())
      expect_that(dim(gmcm$z), equals(c(n, d)))
      # Check K
      expect_that(is.numeric(gmcm$K), is_true())
      expect_that(length(gmcm$K), equals(n))
      # Check theta
      expect_that(is.theta(gmcm$theta), is_true())
    })

  }
}
