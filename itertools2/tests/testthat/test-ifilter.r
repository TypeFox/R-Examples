context("ifilter iterator")

test_that("Apply ifilter to an integer sequence", {
  is_even <- function(x) {
    x %% 2 == 0
  }
  it <- ifilter(is_even, 1:10)

  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 6)
  expect_equal(nextElem(it), 8)
  expect_equal(nextElem(it), 10)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply ifilter to an integer sequence using anonymous function", {
  it <- ifilter(function(x) x %% 2 == 1, 1:10)

  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 5)
  expect_equal(nextElem(it), 7)
  expect_equal(nextElem(it), 9)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply ifilter to a character vector", {
  is_vowel <- function(x) {
    x %in% c('a', 'e', 'i', 'o', 'u')
  }
  it <- ifilter(is_vowel, letters)

  expect_equal(nextElem(it), 'a')
  expect_equal(nextElem(it), 'e')
  expect_equal(nextElem(it), 'i')
  expect_equal(nextElem(it), 'o')
  expect_equal(nextElem(it), 'u')
  expect_error(nextElem(it), "StopIteration")
})
