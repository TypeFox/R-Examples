context("gdiff")

test_that("simple cases are identical with diff", {
  x <- 1:10
  expect_equal(gdiff(x), diff(x))
  expect_equal(gdiff(x, lag = 2), diff(x, lag = 2))
  expect_equal(gdiff(x, differences = 3),
               diff(x, differences = 3))
  expect_equal(gdiff(x, lag = 4, differences = 5),
               diff(x, lag = 4, differences = 5))
})

test_that("can diff a factor", {
  x <- factor(letters)
  expect_equal(gdiff(x, FUN = `!=`),
               rep(TRUE, length(x) - 1))
})
