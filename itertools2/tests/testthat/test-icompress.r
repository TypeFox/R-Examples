context("icompress iterator")

test_that("Apply icompress to an integer sequence", {
  n <- 10
  selectors <- rep(c(F, T), n)
  it <- icompress(1:10, selectors)

  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 6)
  expect_equal(nextElem(it), 8)
  expect_equal(nextElem(it), 10)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply icompress to an integer sequence using anonymous function", {
  n <- 10
  it <- icompress(1:10, rep(c(T, F), n))

  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 5)
  expect_equal(nextElem(it), 7)
  expect_equal(nextElem(it), 9)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply icompress to a character vector", {
  it <- icompress(letters, letters %in% c('a', 'e', 'i', 'o', 'u'))

  expect_equal(nextElem(it), 'a')
  expect_equal(nextElem(it), 'e')
  expect_equal(nextElem(it), 'i')
  expect_equal(nextElem(it), 'o')
  expect_equal(nextElem(it), 'u')
  expect_error(nextElem(it), "StopIteration")
})
