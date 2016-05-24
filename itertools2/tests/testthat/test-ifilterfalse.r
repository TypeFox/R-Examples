context("ifilterfalse iterator")

test_that("Apply ifilterfalse to an integer sequence", {
  is_even <- function(x) {
    x %% 2 == 0
  }
  it <- ifilterfalse(is_even, 1:10)

  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 5)
  expect_equal(nextElem(it), 7)
  expect_equal(nextElem(it), 9)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply ifilterfalse to an integer sequence using anonymous function", {
  it <- ifilterfalse(function(x) x %% 2 == 1, 1:10)

  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 6)
  expect_equal(nextElem(it), 8)
  expect_equal(nextElem(it), 10)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply ifilterfalse to a character vector", {
  is_vowel <- function(x) {
    x %in% c('a', 'e', 'i', 'o', 'u')
  }
  it <- ifilterfalse(is_vowel, letters)

  expect_equal(nextElem(it), 'b')
  expect_equal(nextElem(it), 'c')
  expect_equal(nextElem(it), 'd')
  expect_equal(nextElem(it), 'f')
  expect_equal(nextElem(it), 'g')
  expect_equal(nextElem(it), 'h')
  expect_equal(nextElem(it), 'j')
  expect_equal(nextElem(it), 'k')
  expect_equal(nextElem(it), 'l')
  expect_equal(nextElem(it), 'm')
  expect_equal(nextElem(it), 'n')
  expect_equal(nextElem(it), 'p')
  expect_equal(nextElem(it), 'q')
  expect_equal(nextElem(it), 'r')
  expect_equal(nextElem(it), 's')
  expect_equal(nextElem(it), 't')
  expect_equal(nextElem(it), 'v')
  expect_equal(nextElem(it), 'w')
  expect_equal(nextElem(it), 'x')
  expect_equal(nextElem(it), 'y')
  expect_equal(nextElem(it), 'z')
  expect_error(nextElem(it), "StopIteration")
})
