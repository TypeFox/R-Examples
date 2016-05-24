context("blank-correction")

test_that("grouped blank correction works", {

  file <- system.file("extdata/cary/scans_day_1/", package = "eemR")
  eems <- eem_read(file)

  eems <- eem_remove_blank(eems)

  # These values have been calculated by "hand" in Excel.
  expect_equal(sum(eems[[1]]$x), 14020.78844, tolerance = 0.00001)
  expect_equal(mean(eems[[1]]$x), 1.6038, tolerance = 0.0001)

})


test_that("single blank correction works", {

  file <- system.file("extdata/cary/scans_day_1/sample1.csv", package = "eemR")
  eems <- eem_read(file)

  file <- system.file("extdata/cary/scans_day_1/nano.csv", package = "eemR")
  blank <- eem_read(file)

  eems <- eem_remove_blank(eems, blank)

  # These values have been calculated by "hand" in Excel.
  expect_equal(sum(eems[[1]]$x), 14020.78844, tolerance = 0.00001)
  expect_equal(mean(eems[[1]]$x), 1.6038, tolerance = 0.0001)

})
