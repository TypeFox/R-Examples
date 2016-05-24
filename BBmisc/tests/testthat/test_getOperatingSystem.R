context("getOperatingSystem")

test_that("getOperatingSystem", {
  x = getOperatingSystem()
  expect_true(is.character(x) && length(x) == 1 && nchar(x) > 0)
  x = isWindows()
  expect_true(is.logical(x) && length(x) == 1 && !is.na(x))
  x = isUnix()
  expect_true(is.logical(x) && length(x) == 1 && !is.na(x))
})