context("Testing EMD")

set.seed(1)

test_that("bogus arguments throw error",{
  x <- rnorm(10)
  expect_error(emd(x, noise_strength = 1))
  expect_error(emd(x, ensemble_size = 100))
})

test_that("output of EMD is of correct size and form",{
  x <- rnorm(64)
  imfs <- emd(x, num_imfs = 3)
  expect_identical(dim(imfs), c(64L, 3L))
  expect_identical(class(imfs), c("mts", "ts", "matrix"))
  x <- ts(x, start = 2000, frequency = 12)
  imfs <- emd(x, num_imfs = 3)
  expect_identical(tsp(imfs), tsp(x))
})

test_that("EMD is EEMD",{
  x <- rnorm(10)
  expect_identical(emd(x), eemd(x, ensemble_size = 1, noise_strength = 0))
  x <- rnorm(100)
  expect_identical(emd(x), eemd(x, ensemble_size = 1, noise_strength = 0))
})

test_that("EMD is EEMD",{
  x <- 0
  expect_identical(emd(x), eemd(x, ensemble_size = 1, noise_strength = 0))
  x <- rnorm(100)
  expect_identical(emd(x), eemd(x, ensemble_size = 1, noise_strength = 0))
})

test_that("EMD returns the original series for short series",{
  for (i in 1:3) {
    x <- rnorm(i)
    imfs <- emd(x)
    expect_identical(ncol(imfs), 1L)
    expect_identical(x,c(imfs))
  }
})

test_that("EMD returns the same value each time",{
  for (i in 1:10) {
    x <- rnorm(64)
    expect_identical(emd(x, num_siftings = i), emd(x, num_siftings = i))
    expect_identical(emd(x, S_number = i), emd(x, S_number = i))
  }
})

test_that("sum of imfs equals to original series",{
  x <- rnorm(64)
  expect_equal(rowSums(emd(x)), x)
})

test_that("subsets of IMFs are identical for different num_imfs",{
  x <- rnorm(64)
  imfs3 <- emd(x, num_imfs = 3)
  imfs4 <- emd(x, num_imfs = 4)
  expect_identical(imfs3[, 1:2], imfs4[, 1:2])
})

test_that("num_imfs = 1 returns residual which equals data",{
  x <- rnorm(64)
  imfs <- emd(x, num_imfs = 1)
  expect_identical(c(imfs), x)
})
