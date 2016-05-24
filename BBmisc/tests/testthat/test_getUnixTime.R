context("getUnixTime")

test_that("getUnixTime", {
  x = getUnixTime()
  expect_true(is.integer(x) && length(x) == 1 && !is.na(x))
})
