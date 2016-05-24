library(functools)
context("Lapply()")

x <- data.frame(matrix(runif(10000), ncol = 5))
test_that("Produces the correct output.", {
  expect_equal(Lapply(x, Identity), lapply(x, Identity))
  expect_equal(Lapply(x, sum), lapply(x, sum))
  expect_equal(Lapply(x, mean, trim = 0.2), lapply(x, mean, trim = 0.2))
})

test_that("Produces the correct output type.", {
  expect_is(Lapply(x, Identity), "list")
  expect_is(Lapply(x, sum), "list")
  expect_is(Lapply(x, mean, trim = 0.2), "list")
})

test_that("Produces the correct errors.", {
  expect_error(Lapply(Identity, x))
  expect_error(Lapply(sum, x))
  expect_error(Lapply(mean, x, trim = 0.2))
})
