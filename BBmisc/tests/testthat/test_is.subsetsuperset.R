context("isSubset (and isSuperset)")

test_that("isSubset/isSuperset", {
  x = 1:10
  y = 1:11
  expect_true(isSubset(x, y))
  expect_false(isSubset(y, x))
  expect_true(isSubset(x, y, strict = TRUE))

  x = y
  expect_true(isSubset(x, y))
  expect_false(isSubset(x, y, strict = TRUE))
  expect_true(isSubset(y, x))
  expect_false(isSubset(y, x, strict = TRUE))
})
