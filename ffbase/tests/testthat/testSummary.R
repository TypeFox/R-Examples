library(testthat)

context("summary")

test_that("all ff works",{
  x <- c(TRUE)
  fx <- ff(x)
  expect_equal(all(x), all(fx))

  x <- c(FALSE)
  fx <- ff(x)
  expect_equal(all(x), all(fx))

  x <- c(TRUE, FALSE)
  fx <- ff(x)
  expect_equal(all(x), all(fx))
  
  x <- c(TRUE,NA)
  fx <- ff(x)
  expect_equal(all(x, na.rm=TRUE), all(fx, na.rm=TRUE))
  expect_equal(all(x), all(fx))

  x <- c(FALSE,NA)
  fx <- ff(x)
  expect_equal(all(x, na.rm=TRUE), all(fx, na.rm=TRUE))
  expect_equal(all(x), all(fx))

  x <- c(TRUE, FALSE,NA)
  fx <- ff(x)
  expect_equal(all(x, na.rm=TRUE), all(fx, na.rm=TRUE))
  expect_equal(all(x), all(fx))  
})

test_that("any ff works",{
  x <- c(TRUE)
  fx <- ff(x)
  expect_equal(any(x), any(fx))
  
  x <- c(FALSE)
  fx <- ff(x)
  expect_equal(any(x), any(fx))

  x <- c(TRUE, FALSE)
  fx <- ff(x)
  expect_equal(any(x), any(fx))
  
  x <- c(TRUE, NA)
  fx <- ff(x)
  expect_equal(any(x, na.rm=TRUE), any(fx, na.rm=TRUE))
  expect_equal(any(x), any(fx))
  
  x <- c(FALSE, NA)
  fx <- ff(x)
  expect_equal(any(x, na.rm=TRUE), any(fx, na.rm=TRUE))
  expect_equal(any(x), any(fx))

  x <- c(TRUE, FALSE, NA)
  fx <- ff(x)
  expect_equal(any(x, na.rm=TRUE), any(fx, na.rm=TRUE))
  expect_equal(any(x), any(fx))
})


test_that("Min ff works",{
  x <- runif(100) 
  fx <- ff(x)
  expect_equal(min(x), min(fx))
  
  is.na(x) <- sample(100, 10)
  fx <- ff(x)
  expect_equal(min(x), min(fx))
  expect_equal(min(x, na.rm=TRUE), min(fx, na.rm=TRUE))
})
