context("Check qgmm.marginal")

for (d in seq_len(3) + 1) {
  for (n in trunc(seq(50, 50000, length.out = 20))) {

    sim <- SimulateGMCMData(n = n, theta = rtheta(d = d))
    dimnames(sim$z) <- list(paste0("f", 1:n), paste0("d", 1:d))

    ans <- GMCM:::qgmm.marginal(u = sim$u, theta = rtheta(d = d))

    test_that("qgmm.marginal returns proper format", {
      expect_that(is.numeric(ans), is_true())
      expect_that(dim(ans), equals(c(n, d)))
      expect_that(any(is.na(ans)), is_false())
    })

  }
}
