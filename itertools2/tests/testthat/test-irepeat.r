context("irepeat iterator")

test_that("Indefinite irepeat of integer", {
  it <- irepeat(42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)
})

test_that("irepeat repeats integer fixed number of times", {
  it <- irepeat(42, times=4)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 42)

  expect_error(nextElem(it), "StopIteration")
})

test_that("irepeat repeats data.frame fixed number of times", {
  it <- irepeat(iris, times=5)

  i <- nextElem(it)
  expect_equal(i, iris)

  i <- nextElem(it)
  expect_equal(i, iris)

  i <- nextElem(it)
  expect_equal(i, iris)

  i <- nextElem(it)
  expect_equal(i, iris)

  i <- nextElem(it)
  expect_equal(i, iris)

  expect_error(nextElem(it), "StopIteration")
})
