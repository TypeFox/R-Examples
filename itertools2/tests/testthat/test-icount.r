context("icount iterator: Infinite sequence")

test_that("icount's default values generate the sequence 0, 1, 2, ..", {
  it <- icount()

  i <- nextElem(it)
  expect_equal(i, 0)

  i <- nextElem(it)
  expect_equal(i, 1)

  i <- nextElem(it)
  expect_equal(i, 2)

  i <- nextElem(it)
  expect_equal(i, 3)
})

test_that("icount works with a given initial value", {
  it <- icount(start=42)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 43)

  i <- nextElem(it)
  expect_equal(i, 44)

  i <- nextElem(it)
  expect_equal(i, 45)
})

test_that("icount works with a decimal step size", {
  it <- icount(start=42, step=1.5)

  i <- nextElem(it)
  expect_equal(i, 42)

  i <- nextElem(it)
  expect_equal(i, 43.5)

  i <- nextElem(it)
  expect_equal(i, 45)

  i <- nextElem(it)
  expect_equal(i, 46.5)
})
