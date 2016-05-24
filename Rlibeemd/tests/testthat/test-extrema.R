context("Testing extrema")

set.seed(1)

test_that("bogus arguments throw error",{
  expect_error(extrema("abc"))
  expect_error(extrema(NULL))
})

test_that("length 1 works",{
  ex <- extrema(2)
  expect_identical(c(ex$minima), c(0, 2))
  expect_identical(c(ex$maxima), c(0, 2))
})

test_that("length 2 works",{
  ex <- extrema(c(2, 5))
  expect_identical(c(ex$minima), c(0, 1, 2, 5))
  expect_identical(c(ex$maxima), c(0, 1, 2, 5))
})

test_that("length 3 works",{
  ex <- extrema(c(2, -1, 5))
  expect_identical(c(ex$minima), c(0, 1, 2, 2, -1, 5))
  expect_identical(c(ex$maxima), c(0, 2, 2, 5))
})


test_that("length 3 ts object works",{
  ex <- extrema(ts(c(2, -1, 5), start = 10, frequency = 4))
  expect_identical(c(ex$minima), c(10, 10.25, 10.5, 2, -1, 5))
  expect_identical(c(ex$maxima), c(10, 10.5, 2, 5))
})

test_that("zero series works",{
  ex <- extrema(rep(0,10))
  expect_identical(c(ex$minima), c(0, 9, 0, 0))
  expect_identical(c(ex$maxima), c(0, 9, 0, 0))
})


test_that("all extremas are found correctly",{
  x <- rnorm(128)
  ex <- extrema(x)
  expect_identical(ex$minima[1], 0)
  expect_identical(ex$minima[nrow(ex$minima)], 127)
  expect_identical(ex$maxima[1], 0)
  expect_identical(ex$maxima[nrow(ex$maxima)], 127)
  expect_identical(ex$minima[2:(nrow(ex$minima) - 1),2], 
                    x[1 + ex$minima[2:(nrow(ex$minima) - 1), 1]])
  expect_identical(ex$maxima[2:(nrow(ex$maxima) - 1),2], 
                   x[1 + ex$maxima[2:(nrow(ex$maxima) - 1), 1]])
  
  for (i in 2:(nrow(ex$minima) - 1)) {
    j <- (1 + ex$minima[i, 1])
    expect_true(ex$minima[i,2] <= x[j - 1])
    expect_true(ex$minima[i,2] <= x[j + 1])
  }
  for (i in 2:(nrow(ex$maxima) - 1)) {
    j <- (1 + ex$maxima[i, 1])
    expect_true(ex$maxima[i,2] >= x[j - 1])
    expect_true(ex$maxima[i,2] >= x[j + 1])
  }
})


test_that("check extremas for ts object",{
  x <- ts(rnorm(120), start = 2000, frequency = 4)
  ex <- extrema(x)
  expect_identical(ex$minima[1], 2000)
  expect_identical(ex$minima[nrow(ex$minima)], 2029.75)
  expect_identical(ex$maxima[1], 2000)
  expect_identical(ex$maxima[nrow(ex$maxima)], 2029.75)
  expect_identical(ex$minima[2:(nrow(ex$minima) - 1),2], 
                   x[time(x) %in% (ex$minima[2:(nrow(ex$minima) - 1), 1])])
  expect_identical(ex$maxima[2:(nrow(ex$maxima) - 1),2], 
                   x[time(x) %in% (ex$maxima[2:(nrow(ex$maxima) - 1), 1])])
})
