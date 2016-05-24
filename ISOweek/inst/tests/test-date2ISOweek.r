context("date2ISOweek")

test_that("date2ISOweek returns proper format", {
  expect_that(date2ISOweek("1999-12-31"), is_a("character"))
  expect_that(date2ISOweek("1999-12-31"), matches("[0-9]{4}-W[0-9]{2}-[0-9]{1}"))
})
test_that("date2ISOweek returns proper values", {
  input_chr <- paste(1999:2011, "-12-31", sep="")
  input_date <- as.Date(input_chr)
  result_chr <- c("1999-W52-5", "2000-W52-7", "2002-W01-1", "2003-W01-2", "2004-W01-3", 
                  "2004-W53-5", "2005-W52-6", "2006-W52-7", "2008-W01-1", "2009-W01-3",
                  "2009-W53-4", "2010-W52-5", "2011-W52-6")
  expect_that(date2ISOweek(input_date), equals(result_chr))
  expect_that(date2ISOweek(input_chr), equals(result_chr))
})
test_that("date2ISOweek handles NAs", {
  expect_that(date2ISOweek(NA), equals(NA_character_))
  expect_that(date2ISOweek(NA_character_), equals(NA_character_))
  input_chr <- paste(1999:2011, "-12-31", sep="")
  result_chr <- c("1999-W52-5", "2000-W52-7", "2002-W01-1", "2003-W01-2", "2004-W01-3", 
                  "2004-W53-5", "2005-W52-6", "2006-W52-7", "2008-W01-1", "2009-W01-3",
                  "2009-W53-4", "2010-W52-5", "2011-W52-6")
  idx_NA <- seq(from =2, to = length(input_chr), by = 2)
  input_chr_NA <- input_chr
  input_chr_NA[idx_NA] <- NA_character_
  input_date_NA <- as.Date(input_chr_NA)
  result_chr_NA <- result_chr
  result_chr_NA[idx_NA] <- NA
  expect_that(date2ISOweek(input_date_NA), equals(result_chr_NA))
  expect_that(date2ISOweek(input_chr_NA), equals(result_chr_NA))
})
test_that("date2ISOweek handles zero length vectors", {
  input_date <- as.Date("1999-12-31")
  expect_that(length(date2ISOweek(input_date[0])), equals(0))
})
test_that("date2ISOweek stops on invalid parameters", {
  expect_that(date2ISOweek(31), throws_error())
  expect_that(date2ISOweek(31.1), throws_error())
  expect_that(date2ISOweek(FALSE), throws_error())
  expect_that(date2ISOweek("31"), throws_error())
  expect_that(date2ISOweek("12/31/1999"), throws_error())
})
