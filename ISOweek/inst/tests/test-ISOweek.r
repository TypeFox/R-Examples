context("ISOweek")

test_that("ISOweek returns proper format", {
  expect_that(ISOweek("1999-12-31"), is_a("character"))
  expect_that(ISOweek("1999-12-31"), matches("[0-9]{4}-W[0-9]{2}"))
})
test_that("ISOweek returns proper values", {
  input_chr <- paste(1999:2011, "-12-31", sep="")
  input_date <- as.Date(input_chr)
  result_chr <- c("1999-W52", "2000-W52", "2002-W01", "2003-W01", "2004-W01", 
                  "2004-W53", "2005-W52", "2006-W52", "2008-W01", "2009-W01",
                  "2009-W53", "2010-W52", "2011-W52")
  expect_that(ISOweek(input_date), equals(result_chr))
  expect_that(ISOweek(input_chr), equals(result_chr))
})
test_that("ISOweek handles NAs", {
  expect_that(ISOweek(NA), equals(NA_character_))
  expect_that(ISOweek(NA_character_), equals(NA_character_))
  input_chr <- paste(1999:2011, "-12-31", sep="")
  result_chr <- c("1999-W52", "2000-W52", "2002-W01", "2003-W01", "2004-W01", 
                  "2004-W53", "2005-W52", "2006-W52", "2008-W01", "2009-W01",
                  "2009-W53", "2010-W52", "2011-W52")
  idx_NA <- seq(from =2, to = length(input_chr), by = 2)
  input_chr_NA <- input_chr
  input_chr_NA[idx_NA] <- NA_character_
  input_date_NA <- as.Date(input_chr_NA)
  result_chr_NA <- result_chr
  result_chr_NA[idx_NA] <- NA
  expect_that(ISOweek(input_date_NA), equals(result_chr_NA))
  expect_that(ISOweek(input_chr_NA), equals(result_chr_NA))
})
test_that("ISOweek handles zero length vectors", {
  input_date <- as.Date("1999-12-31")
  expect_that(length(ISOweek(input_date[0])), equals(0))
})
test_that("ISOweek stops on invalid parameters", {
  expect_that(ISOweek(31), throws_error())
  expect_that(ISOweek(31.1), throws_error())
  expect_that(ISOweek(FALSE), throws_error())
  expect_that(ISOweek("31"), throws_error())
  expect_that(ISOweek("12/31/1999"), throws_error())
})
