context("dframe")

df <- sample_data_1

test_that("dframe basic use without lat/long vars works", {
  skip_on_cran()

  aa <- suppressMessages(dframe(df))

  expect_is(aa, "data.frame")
  expect_is(aa, "dframe")
})

test_that("dframe fails well", {
  skip_on_cran()

  expect_error(dframe(),
               "argument \"x\" is missing")
  expect_error(dframe("things"),
               "no 'dframe' method for")
})
