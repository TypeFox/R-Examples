context("Test spread_to_upfront")

## Test. In the following test cases we compare the results of 
## our upront functions for different data with results from markit.com using
## the same data

test_that("test for spread_to_upfront", {
  
  ## test case to see if our function gives the same result as markit.com
  ## all cases use data from Xerox Corporation for 2014-04-22. 
  
  load("test_spread_to_upfront.RData")
  
  ## actual upfront value from markit.com for Xerox Corporation for 2014-04-22.
  
  ## truth.1 <- 18624
  
  result.1 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-22"),
                                           currency = "USD",
                                           tenor    = 5,
                                           spread   = 105.8,
                                           coupon   = 100,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE),
                                isPriceClean = FALSE)
  
  ## Note: test case passes when rounded to the nearest tenth.
  ## Difference of 3.39 (0.02094072 %) from actual value
  
  expect_that(round(result.1), equals(round(truth.1)))
  
  ## In the following test cases we want to check if the result of changing
  ## certain variables is the same as the results from markit.com results.
  
  ## test case where spread is equal to the coupon 
  
  ## markit.com value
  ## truth.2 <- -9444
  ## calculated value
  
  result.2 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-22"),
                                           currency = "USD",
                                           tenor    = 5,
                                           spread   = 100,
                                           coupon   = 100,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  ## comparing the results with markit data
  ## Note: test case passes when the values are rounded off till the nearest whole
  ## number
  ## difference of 0.444 from actual number
  
  expect_that(round(result.2), equals(truth.2))
  
  
  ## Effect on upfront of an increase in coupon rate (by 100 basis points). 
  
  ## actual value
  ## truth.3 <- -474755
  ## calculated value
  
  result.3 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-22"),
                                           currency = "USD",
                                           tenor    = 5,
                                           spread   = 105.8,
                                           coupon   = 200,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  ## comparing the results with markit data
  ## Note: test case passes when results are rounded off the nearest 100
  ## difference of $14 (0.00294902 %) from actual number
  
  expect_that(round(result.3), equals(round(truth.3)))
  
  
  ## Effect on upfront of an decrease in coupon rate (by 50 basis points). 
  
  ## actual value
  ## truth.4 <- 265313
  ## calculated value
  
  result.4 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-22"),
                                           currency = "USD",
                                           tenor    = 5,
                                           spread   = 105.8,
                                           coupon   = 50,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  ## comparing the results with markit data
  ## Note: test case passes when results are rounded off the nearest 100
  ## difference of $8.408048 (0.003169249 %) from actual number
  
  expect_that(round(result.4), equals(round(truth.4)))
  
  
  ## Effect on upfront of a decrease in trade date (by one week)
  
  ## actual value
  ## truth.5 <- 20718
  ## calculated value
  
  result.5 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-15"),
                                           currency = "USD",
                                           maturity = as.Date("2019-06-20"),
                                           spread   = 105.8,
                                           coupon   = 100,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  ## comparing the results with markit data
  ## Note: test case passes when rounded to the nearest 10
  ## difference of $0.9439034 (0.004556178 %) from actual number
  
  expect_that(round(result.5), equals(round(truth.5)))
  
  
  ## Effect on upfront of an increase in the trade date (by one week)
  
  ## actual value
  ## truth.6 <- 16582
  #calculated value
  
  result.6 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-29"),
                                           currency = "USD",
                                           maturity = as.Date("2019-06-20"),
                                           spread   = 105.8,
                                           coupon   = 100,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  ## comparing the results with markit data
  ## Note: test case passes when rounded of to the nearest tenth
  
  expect_that(round(result.6), equals(round(truth.6)))
  
  
  ## Effect on upfront of a decrease in the maturity date (by one quarter)
  
  ## actual value when maturity is 2019-09-20 instead of 2019-06-20
  ## truth.7 <- 17395
  #calculated value
  
  result.7 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-22"),
                                           currency = "USD",
                                           maturity = as.Date("2019-03-20"),
                                           spread   = 105.8,
                                           coupon   = 100,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  ## comparing the results with markit data
  
  expect_that(round(result.7), equals(round(truth.7)))
  
  
  ## Effect on upfront of an increase in the maturity date (by one quarter)
  
  ## actual value
  ## truth.8 <- 19836
  ## calculated value
  
  result.8 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-22"),
                                           currency = "USD",
                                           maturity = as.Date("2019-09-20"),
                                           spread   = 105.8,
                                           coupon   = 100,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  
  ## comparing the results with markit data
  ## Note: test case passes when rounded off to the nearest 10.
  
  expect_that(round(result.8), equals(round(truth.8)))
  
  ## Effect on upfront of an increase in spread (by 50 basis points)
  
  ## actual value
  ## truth.9 <- 254985
  ## calculated value
  
  result.9 <- spread_to_upfront(data.frame(date     = as.Date("2014-04-22"),
                                           currency = "USD",
                                           tenor    = 5,
                                           spread   = 155.8,
                                           coupon   = 100,
                                           recovery = 0.4,
                                           stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  ## comparing the results with markit data
  ## Note: test case passes when rounded to nearest 1000
  
  expect_that(round(result.9), equals(round(truth.9)))
  
  
  ## Effect on upfront of an decrease in spread (by 50 basis points)
  
  ## actual value
  ## truth.10 <- -227912
  ## calculated value
  
  result.10 <- spread_to_upfront(data.frame(date = as.Date("2014-04-22"),
                                            currency = "USD",
                                            tenor = 5,
                                            spread = 55.8,
                                            coupon = 100,
                                            recovery = 0.4,
                                            stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  
  ## comparing the results with markit data
  ## Note: test case passes when rounded off to the nearest 1000
  
  expect_that(round(result.10), equals(round(truth.10)))
  
  
  ## Effect on upfront when trade date = maturity date (September 20, 2013)
  
  ## actual value
  ## truth.11 <- 0
  ## calculated value
  ## 
  
  result.11 <- spread_to_upfront(data.frame(date = as.Date("2013-09-20"),
                                            currency = "USD",
                                            maturity = as.Date("2013-09-20"),
                                            spread = 105.8,
                                            coupon = 100,
                                            recovery = 0.4,
                                            stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  
  # comparing the results with markit data
  ## Note: test case passes when rounded off to the nearest 1000
  ## expect_that(round(result.11), equals(round(truth.11)))
  
  
  ## save(truth.1, truth.2, truth.3, truth.4, truth.5, truth.6, truth.7, 
  ##     truth.8, truth.9, truth.10, truth.11, file="upfront.test.RData")
  
  ## test case to see an the effect on upfront when trade date = roll over date
  
  ## test case to see upfront payment when trade date is one day away from 
  ## maturity date
  
  ## test for different Japanese dates
  
  result.13 <- spread_to_upfront(data.frame(date = as.Date("2009-03-18"),
                                            currency = "JPY",
                                            maturity = as.Date("2014-03-20"),
                                            spread = 105.8,
                                            coupon = 100,
                                            recovery = 0.35,
                                            stringsAsFactors = FALSE), isPriceClean = FALSE)
  
  truth.13 <- 3487
  
  expect_true(abs(result.13 - truth.13) < 5)
  
  x.1 <- data.frame(date = as.Date(c("2014-04-15", "2014-04-22")),  
                    tenor = c(5, 5), 
                    coupon = c(500, 100), 
                    spread = c(2785.8889, 99),
                    currency = c("EUR", "EUR"),
                    recovery = c(0.4, 0.4),
                    stringsAsFactors = FALSE)
  
  x.2 <- data.frame(date = as.Date(c("2014-04-15", "2014-04-22")),
                    maturity = as.Date(c("2019-06-20", "2019-06-20")),
                    coupon = c(500, 100), 
                    spread = c(2785.8889, 99),
                    currency = c("EUR", "EUR"),
                    recovery = c(0.4, 0.4),
                    stringsAsFactors = FALSE)
  
  result.1 <- spread_to_upfront(x = x.1, tenor.var = "tenor")
  expect_that(round(result.1, -2), equals(c(4412500, round(-14368, -2))))
  
  result.2 <- spread_to_upfront(x.2, tenor = "tenor")
  expect_that(round(result.2, -2), equals(c(4412500, round(-14368, -2))))
})