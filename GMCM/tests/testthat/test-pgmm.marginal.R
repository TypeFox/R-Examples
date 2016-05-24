context("Check pgmm.marginal")

for (d in seq_len(3) + 1) {
  for (n in trunc(seq(50, 50000, length.out = 20))) {

    sim <- SimulateGMMData(n = n, theta = rtheta(d = d))
    dimnames(sim$z) <- list(paste0("f", 1:n), paste0("d", 1:d))

    ans <- GMCM:::pgmm.marginal(z = sim$z, theta = rtheta(d = d))

    test_that("pgmm.marginal returns proper format", {
      expect_that(is.numeric(ans), is_true())
      expect_that(dim(ans), equals(c(n, d)))
      expect_that(min(ans), is_more_than(0 - .Machine$double.eps))
      expect_that(max(ans), is_less_than(1 + .Machine$double.eps))
    })

  }
}
