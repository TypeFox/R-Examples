context("iproduct iterator")

test_that("iproduct constructs the Cartesian product of two unnamed numeric vectors", {
  it <- iproduct(1:3, 4:6)
  expect_equal(nextElem(it), list(1, 4))
  expect_equal(nextElem(it), list(1, 5))
  expect_equal(nextElem(it), list(1, 6))
  expect_equal(nextElem(it), list(2, 4))
  expect_equal(nextElem(it), list(2, 5))
  expect_equal(nextElem(it), list(2, 6))
  expect_equal(nextElem(it), list(3, 4))
  expect_equal(nextElem(it), list(3, 5))
  expect_equal(nextElem(it), list(3, 6))
  
  expect_error(nextElem(it), "StopIteration")
})

test_that("iproduct constructs the Cartesian product of three unnamed numeric vectors", {
  it <- iproduct(1:2, 3:4, 5:6)
  expect_equal(nextElem(it), list(1, 3, 5))
  expect_equal(nextElem(it), list(1, 3, 6))
  expect_equal(nextElem(it), list(1, 4, 5))
  expect_equal(nextElem(it), list(1, 4, 6))
  expect_equal(nextElem(it), list(2, 3, 5))
  expect_equal(nextElem(it), list(2, 3, 6))
  expect_equal(nextElem(it), list(2, 4, 5))
  expect_equal(nextElem(it), list(2, 4, 6))

  expect_error(nextElem(it), "StopIteration")
})

test_that("iproduct constructs the Cartesian product of two named numeric vectors", {
  it <- iproduct(a=1:3, b=4:6)
  expect_equal(nextElem(it), list(a=1, b=4))
  expect_equal(nextElem(it), list(a=1, b=5))
  expect_equal(nextElem(it), list(a=1, b=6))
  expect_equal(nextElem(it), list(a=2, b=4))
  expect_equal(nextElem(it), list(a=2, b=5))
  expect_equal(nextElem(it), list(a=2, b=6))
  expect_equal(nextElem(it), list(a=3, b=4))
  expect_equal(nextElem(it), list(a=3, b=5))
  expect_equal(nextElem(it), list(a=3, b=6))
  
  expect_error(nextElem(it), "StopIteration")
})

test_that("iproduct constructs the Cartesian product of three named numeric vectors", {
  it <- iproduct(a=1:2, b=3:4, c=5:6)
  expect_equal(nextElem(it), list(a=1, b=3, c=5))
  expect_equal(nextElem(it), list(a=1, b=3, c=6))
  expect_equal(nextElem(it), list(a=1, b=4, c=5))
  expect_equal(nextElem(it), list(a=1, b=4, c=6))
  expect_equal(nextElem(it), list(a=2, b=3, c=5))
  expect_equal(nextElem(it), list(a=2, b=3, c=6))
  expect_equal(nextElem(it), list(a=2, b=4, c=5))
  expect_equal(nextElem(it), list(a=2, b=4, c=6))

  expect_error(nextElem(it), "StopIteration")
})

