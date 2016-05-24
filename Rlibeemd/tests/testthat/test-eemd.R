context("Testing EEMD")

set.seed(1)

test_that("bogus arguments throw error",{
  expect_error(eemd("abc"))
  expect_error(eemd(1:3, noise_strength = -1, threads = 1))
  expect_error(eemd(1:3, num_siftings = 0, S_number = 0, threads = 1))
  expect_error(eemd(1:3, num_imfs = -1, threads = 1))
  expect_error(eemd(1:3, num_siftings = -3, threads = 1))
  expect_error(eemd(1:3, ensemble_size = 0, threads = 1))
  expect_error(eemd(1:3, ensemble_size = 10, noise_strength = 0, threads = 1))
  expect_error(eemd(1:3, ensemble_size = 1, noise_strength = 1, threads = 1))
  expect_error(eemd(1:3, num_imfs = "lots", threads = 1))
})

test_that("series full of zeroes should produce only zeroes",{
  x <- ts(numeric(64))
  imfs <- eemd(x, ensemble_size = 10, threads = 1)
  expect_true(all(imfs == 0))
})

test_that("residual is close to the disturbed signal",{
  n <- 100
  x <- seq(1, 10, length.out = n) ^ 2 + rnorm(n, sd = 0.5)
  imfs <- eemd(x, threads = 1)
  residual <- imfs[, 6]
  expect_equal(residual[11:90], x[11:90], tol = 0.1)
})

test_that("different seeds give different results",{
  x <- rnorm(64)
  expect_false(isTRUE(all.equal(eemd(x, rng_seed = 1, threads = 1), 
                                eemd(x, rng_seed = 2, threads = 1))))
})

test_that("identical seeds give equal results",{
  x <- rnorm(64)
  expect_identical(eemd(x, rng_seed = 1, threads = 1), 
               eemd(x, rng_seed = 1, threads = 1))
})

test_that("subsets of IMFs are identical for different num_imfs",{
  x <- rnorm(64)
  imfs3 <- eemd(x, num_imfs = 3, rng_seed = 1, threads = 1)
  imfs4 <- eemd(x, num_imfs = 4, rng_seed = 1, threads = 1)
  expect_equal(imfs3[, 1:2], imfs4[, 1:2])
})
