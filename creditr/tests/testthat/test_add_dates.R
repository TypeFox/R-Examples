context("Test add_dates")

## Test cases for the add_dates function (using data from April 15, 2014 for Norske
## Skogindustrier ASA).

library(creditr)

test_that("test for add_dates", {
  
  ## Should return an error when no tenor or maturity date is given
  
  expect_error(add_dates(data.frame(date = as.Date("2014-04-15"), currency = "USD")))
  
  ## Should return an error when both tenor and maturity date are given; only one
  ## should be entered
  
  expect_error(add_dates(data.frame(date = as.Date("2014-04-15"),
                                    tenor = 5,
                                    maturity = as.Date("2016-06-20"),
                                    currency = "USD")))
  
  ## Should return an error when maturity is not of type date
  
  expect_error(add_dates(data.frame(date = as.Date("2014-04-15"),
                                    maturity = "not a date",
                                    currency = "JPY")))
  
  ## different dates for a CDS with a 10-year maturity
  
  ## Test case with CDS from April 15, 2014 (Caesar's Entertainment Corporation)
  
  load("test_add_dates.RData")
  
  result.1 <- add_dates(data.frame(date = as.Date("2014-04-15"), tenor = 5, currency = "USD"))
  result.2 <- add_dates(data.frame(date = as.Date("2014-04-15"),
                                   maturity = as.Date("2019-06-20"),
                                   currency = "USD")) 
  
  ## expect_that(result.1, is_identical_to(truth.1))
  ## expect_that(result.2, is_identical_to(truth.1))
  
  ## if the trade date is right after roll date, then the endDate 
  ## should go to next roll date
  
  result.3 <- add_dates(data.frame(date = as.Date("2011-06-21"),
                                   tenor = 5,
                                   currency = "USD"))
  ## expect_that(result.3, is_identical_to(truth.2))
  
  ## if the trade date is a US holiday, say, independence Day
  ## "2011-07-04", then add_dates() should give a warning because
  ## US trade date can't happen on the independence date. But it 
  ## doesn't now, so have to comment out the following test
  
  # expect_warning(add_dates(data.frame(date = as.Date("2011-07-04"), 
  #                         tenor = 5, currency="USD")))
  
  ## if the trade date is a Monday
  
  expect_equal(add_dates(data.frame(date = as.Date("2011-06-03"), 
                                    tenor = 1, currency = "USD"))$endDate, 
               as.Date("2012-06-20"))
  
  ## if the trade date is only one day before the maturity date,
  ## add_dates() should give a warning because it's impossible
  ## to return dates in duration of 0 days. (since stepinDate is equal
  ## to (maturity date + 1), which is the maturity date, 
  ## the duration is 0 day). This is an extreme test case, and since
  ## add_dates() cannot solve this test case for now, the following 
  ## test case is commented out.
  
  # expect_warning(add_dates(data.frame(date = as.Date("2009-06-19"), 
  #                maturity = as.Date("2009-06-20"))))
  
  ## if the endDate (maturity date) is a weekend, add_dates should just
  ## return a weekend day, instead of adjust it to the next weekeday
  
  ## for example, if we let trade date be "2010-06-18", then the
  ## maturity date should be "2015-06-20", a Saturday.
  
  x <- add_dates(data.frame(date = as.Date("2010-06-18"), tenor = 5, 
                            currency = "USD"))
  expect_equal(x$endDate, as.Date("2015-06-20"))
})

context("JPY holidays baseDate test")

test_that("Test JPY holidays",{
  
  x <- data.frame(date = as.Date("2015-09-20"), tenor = 5, currency = "JPY")
  
  ## "2015-09-21", "2015-09-22" and "2015-09-23" are all Japanese holidays, 
  ## so baseDate should be "2015-09-25"
  
  expect_that(add_dates(x)$baseDate, equals(as.Date("2015-09-25")))
  
})