context("isProperlyNamed")

test_that("isProperlyNamed", {
  expect_true(isProperlyNamed(list()))
  expect_true(isProperlyNamed(list(x=1)))
  expect_true(isProperlyNamed(list(x=1, y=2)))
  expect_true(!isProperlyNamed(list(1,2)))
  xs = list(1,2)
  names(xs)[1] = "a"
  expect_true(!isProperlyNamed(xs))
  names(xs)[2] = "b"
  expect_true(isProperlyNamed(xs))
})