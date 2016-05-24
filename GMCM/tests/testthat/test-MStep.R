context("Check MStep")

for (n in trunc(seq(10, 1000, by = 100))) {
  for (d in seq_len(4) + 1) {

    x <- SimulateGMMData(n = n, d = d,
                         theta = meta2full(c(0.5, 1, 1, 0.5), d = d))$z
    kappa <- matrix(runif(n*d), n, d)
    kappa <- kappa/rowSums(kappa)

    theta <- GMCM:::MStep(x, kappa, meta.special.case = FALSE)

    test_that("MStep returns proper format", {
      expect_that(is.theta(theta), is_true())
    })

  }
}
