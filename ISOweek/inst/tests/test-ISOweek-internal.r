context("Internal functions")

test_that("weekday0 returns proper values", {
  input_chr <- paste(1999:2011, "-12-31", sep="")
  input_date <- as.Date(input_chr)
  result_int <- as.integer(c(5, 7, 1, 2, 3, 5, 6, 7, 1, 3, 4, 5, 6)) - 1L
  expect_that(weekday0(input_date), equals(result_int))
  expect_that(weekday0(input_chr), equals(result_int))
})

test_that("thursday0 returns proper values", {
  expect_that(thursday0(Sys.Date()), is_a("Date"))
  expect_that(weekday0(thursday0("1999-12-31")), equals(3L))
  expect_that(abs(difftime(thursday0(Sys.Date()), Sys.Date(), units = "days")) <= 3, is_true())
})

test_that("year0 returns proper values", {
  expect_that(year0(Sys.Date()), is_a("integer"))
  expect_that(year0(c("1999-12-31", "2000-01-01")), equals(c(1999, 2000)))
})
