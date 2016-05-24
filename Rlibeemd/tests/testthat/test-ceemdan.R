context("Testing CEEMDAN")

set.seed(1)

test_that("bogus arguments throw error",{
  expect_error(ceemdan("abc"))
  expect_error(ceemdan(1:3, noise_strength = -1, threads = 1))
  expect_error(ceemdan(1:3, num_siftings = 0, S_number = 0, threads = 1))
  expect_error(ceemdan(1:3, num_imfs = -1, threads = 1))
  expect_error(ceemdan(1:3, num_siftings = -3, threads = 1))
  expect_error(ceemdan(1:3, ensemble_size = 0, threads = 1))
  expect_error(ceemdan(1:3, num_imfs = "lots", threads = 1))
})

test_that("series full of zeroes should produce only zeroes",{
  x <- numeric(64)
  imfs <- ceemdan(x, ensemble_size = 10, threads = 1)
  expect_true(all(imfs == 0))
})

test_that("residual is close to the disturbed signal",{
  n <- 100
  x <- seq(1, 10, length.out = n) ^ 2 + rnorm(n, sd = 0.5)
  imfs <- ceemdan(x, threads = 1)
  residual <- imfs[, 6]
  expect_equal(residual[11:90], x[11:90], tol = 0.1)
})

test_that("different seeds give different results",{
  x <- rnorm(64)
  expect_false(isTRUE(all.equal(ceemdan(x, rng_seed = 1, threads = 1), 
                                ceemdan(x, rng_seed = 2, threads = 1))))
})

test_that("identical seeds give equal results",{
  x <- rnorm(64)
  expect_identical(ceemdan(x, rng_seed = 1, threads = 1), 
                   ceemdan(x, rng_seed = 1, threads = 1))
})

test_that("subsets of IMFs are identical for different num_imfs",{
  x <- rnorm(64)
  imfs3 <- ceemdan(x, num_imfs = 3, rng_seed = 1, threads = 1)
  imfs4 <- ceemdan(x, num_imfs = 4, rng_seed = 1, threads = 1)
  expect_equal(imfs3[, 1:2], imfs4[, 1:2])
})


test_that("num_imfs = 1 returns residual which equals data",{
  x <- rnorm(64)
  imfs <- ceemdan(x, num_imfs = 1, threads = 1)
  expect_identical(c(imfs), x)
})


test_that("sum of imfs equals to original series",{
  x <- rnorm(64)
  expect_equal(rowSums(ceemdan(x, threads = 1)), x)
})
