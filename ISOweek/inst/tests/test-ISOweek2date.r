context("ISOweek2date")

test_that("ISOweek2date returns proper format", {
  expect_that(ISOweek2date("1999-W01-1"), is_a("Date"))
})
test_that("ISOweek2date returns proper values", {
  input_chr <- c("1999-W52-5", "2000-W52-7", "2002-W01-1", "2003-W01-2", "2004-W01-3", 
                  "2004-W53-5", "2005-W52-6", "2006-W52-7", "2008-W01-1", "2009-W01-3",
                  "2009-W53-4", "2010-W52-5", "2011-W52-6")
  result_chr <- paste(1999:2011, "-12-31", sep="")
  result_date <- as.Date(result_chr)
  expect_that(ISOweek2date(input_chr), equals(result_date))
})
test_that("ISOweek2date handles NAs", {
  expect_that(ISOweek2date(NA), equals(as.Date(NA)))
  expect_that(ISOweek2date(NA_character_), equals(as.Date(NA)))
  input_chr <- c("1999-W52-5", "2000-W52-7", "2002-W01-1", "2003-W01-2", "2004-W01-3", 
                  "2004-W53-5", "2005-W52-6", "2006-W52-7", "2008-W01-1", "2009-W01-3",
                  "2009-W53-4", "2010-W52-5", "2011-W52-6")
  result_chr <- paste(1999:2011, "-12-31", sep="")
  result_date <- as.Date(result_chr)
  idx_NA <- seq(from =2, to = length(input_chr), by = 2)
  input_chr_NA <- input_chr
  input_chr_NA[idx_NA] <- NA_character_
  result_chr_NA <- result_chr
  result_chr_NA[idx_NA] <- NA
  result_date_NA <- as.Date(result_chr_NA)
  expect_that(ISOweek2date(input_chr_NA), equals(result_date_NA))
})
test_that("ISOweek2date handles zero length vectors", {
  input_chr <- c("1999-W52-5")
  expect_that(length(ISOweek2date(input_chr[0])), equals(0))
})
test_that("ISOweek2date stops on invalid parameters", {
  expect_that(ISOweek2date(31), throws_error())
  expect_that(ISOweek2date(31.1), throws_error())
  expect_that(ISOweek2date(FALSE), throws_error())
  expect_that(ISOweek2date("31"), throws_error())
  expect_that(ISOweek2date("12/31/1999"), throws_error())
  expect_that(ISOweek2date("1999-W00-1"), throws_error())
  expect_that(ISOweek2date("1999-W54-1"), throws_error())
  expect_that(ISOweek2date("1999-W01-0"), throws_error())
  expect_that(ISOweek2date("1999-W01-8"), throws_error())
})
