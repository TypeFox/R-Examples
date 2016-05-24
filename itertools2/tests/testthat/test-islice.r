context("islice iterator")

test_that("islice of integer sequence with no end specified", {
  it <- islice(1:5, start=2)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 5)
  expect_error(nextElem(it), "StopIteration")
})

test_that("islice of integer sequence with start and end specified", {
  it <- islice(1:10, start=2, end=5)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 5)
  expect_error(nextElem(it), "StopIteration")
})

test_that("islice of integer sequence with start, end, and step specified", {
  it <- islice(1:10, start=2, end=9, step=2)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 6)
  expect_equal(nextElem(it), 8)
  expect_error(nextElem(it), "StopIteration")
})


