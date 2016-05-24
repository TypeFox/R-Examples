context("iseq family of iterators")

test_that("iseq with default parameters yields only 1 and then StopIteration", {
  it <- iseq()
  expect_equal(nextElem(it), 1)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq with only from and to specified", {
  it <- iseq(from=2, to=5.5)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 5)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq with by specified", {
  it <- iseq(from=2, to=3.5, by=0.6)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 2.6)
  expect_equal(nextElem(it), 3.2)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq with length_out specified", {
  it <- iseq(from=2, to=3.5, length_out=5)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 2.375)
  expect_equal(nextElem(it), 2.75)
  expect_equal(nextElem(it), 3.125)
  expect_equal(nextElem(it), 3.5)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq with along_with specified", {
  it <- iseq(from=2, to=3.5, along_with=1:5)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 2.375)
  expect_equal(nextElem(it), 2.75)
  expect_equal(nextElem(it), 3.125)
  expect_equal(nextElem(it), 3.5)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq for decreasing sequence with from and to specified", {
  it <- iseq(from=2, to=-3.5)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 0)
  expect_equal(nextElem(it), -1)
  expect_equal(nextElem(it), -2)
  expect_equal(nextElem(it), -3)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq_len generates a finite sequence of integers", {
  it <- iseq_len(4)
  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 4)
  expect_error(nextElem(it), "StopIteration")
})

test_that("First element of iseq_len with length 0 yields StopIteration", {
  it <- iseq_len(0)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq_along's generate a finite sequence of integers from a vector", {
  it <- iseq_along(1:4)
  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 4)
  expect_error(nextElem(it), "StopIteration")
})

test_that("iseq_along's generate a finite sequence of integers from a data.frame", {
  it <- iseq_along(iris)
  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 5)
  expect_error(nextElem(it), "StopIteration")
})

test_that("First element of iseq_along applied to vector of length 0 yields StopIteration", {
  it <- iseq_along(numeric(0))
  expect_error(nextElem(it), "StopIteration")
})
