context("Check estimator output")

test_that("Compare Stein and Stein+ISO to Lin-Perlman example output", {
  set.seed(1)
  p <- 10
  n <- 10
  l <- c(50.91, 39.66, 14.96, 11.88, 8.55, 6.33, 2.40, 1.99, 0.52, 0.01) / n

  phi_iso <- iso_eig(l, n)
  ref_iso <- c(2.1778, 2.1778, 0.8606, 0.8606, 0.8606, 0.8606, 0.5671, 0.5671, 0.4129, 0.0107)
  expect_equal(phi_iso, ref_iso, tolerance=1e-4)

  phi_haff <- haff_eig(l, n)
  ref_haff <- c(2.246, 2.246, 0.9908, 0.9908, 0.9908, 0.9908, 0.7476, 0.7476, 0.4129, 0.0107)
  expect_equal(phi_haff, ref_haff, tolerance=1e-4)

  H <- svd(rWishart(1, n, diag(p))[,,1])$u
  S <- H %*% (l * t(H))

  Sigma_iso <- iso_cov(S, n)
  expect_equal(eigen(Sigma_iso)$val, ref_iso, tolerance=1e-4)

  Sigma_haff <- haff_cov(S, n)
  expect_equal(eigen(Sigma_haff)$val, ref_haff, tolerance=1e-4)
})
