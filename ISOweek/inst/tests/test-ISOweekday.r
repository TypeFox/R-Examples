context("ISOweekday")

test_that("ISOweekday returns integer values", {
  expect_that(ISOweekday("1999-12-31"), is_a("integer"))
})
test_that("ISOweekday returns proper values", {
  input_chr <- paste(1999:2011, "-12-31", sep="")
  input_date <- as.Date(input_chr)
  result_int <- as.integer(c(5, 7, 1, 2, 3, 5, 6, 7, 1, 3, 4, 5, 6))
  expect_that(ISOweekday(input_date), equals(result_int))
  expect_that(ISOweekday(input_chr), equals(result_int))
})
test_that("ISOweekday handles NAs", {
  expect_that(ISOweekday(NA), equals(NA_integer_))
  expect_that(ISOweekday(NA_character_), equals(NA_integer_))
  input_chr <- paste(1999:2011, "-12-31", sep="")
  result_int <- as.integer(c(5, 7, 1, 2, 3, 5, 6, 7, 1, 3, 4, 5, 6))
  idx_NA <- seq(from =2, to = length(input_chr), by = 2)
  input_chr_NA <- input_chr
  input_chr_NA[idx_NA] <- NA_character_
  input_date_NA <- as.Date(input_chr_NA)
  result_int_NA <- result_int
  result_int_NA[idx_NA] <- NA
  expect_that(ISOweekday(input_date_NA), equals(result_int_NA))
  expect_that(ISOweekday(input_chr_NA), equals(result_int_NA))
})
test_that("ISOweekday handles zero length vectors", {
  input_date <- as.Date("1999-12-31")
  expect_that(length(ISOweekday(input_date[0])), equals(0))
})
test_that("ISOweekday stops on invalid parameters", {
  expect_that(ISOweekday(31), throws_error())
  expect_that(ISOweekday(31.1), throws_error())
  expect_that(ISOweekday(FALSE), throws_error())
  expect_that(ISOweekday("31"), throws_error())
  expect_that(ISOweekday("12/31/1999"), throws_error())
})
