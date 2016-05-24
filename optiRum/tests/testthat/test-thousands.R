context("thousands")
library(plyr)
test_that("thousands format always adds k", {
    expect_equal(thousands(1000), "1k")
    expect_equal(thousands(1e+06), "1,000k")
    expect_equal(thousands(1e+09), "1,000,000k")
})

test_that("thousands_format format always adds k", {
  expect_equal(thousands_format()(1000), "1k")
  expect_equal(thousands_format()(1e+06), "1,000k")
  expect_equal(thousands_format()(1e+09), "1,000,000k")
})
