context("Test upfront_to_spread")

## upfront_to_spread() test case with data from Xerox corporation

test_that("test for upfront_to_spread", {
  ## actual spread value from Bloomberg for Xerox Corp.
   truth.1 <- 105.85
  
  load("test_upfront_to_spread.RData")

  result.1 <- upfront_to_spread(x = data.frame(date = as.Date("2014-04-22"),
                                               currency = "USD",
                                               coupon = 100,
                                               recovery = .4,
                                               tenor = 5,
                                               upfront = 18624),
                                isPriceClean = FALSE,
                                notional = 10000000,
                                payAccruedAtStart = FALSE,
                                payAccruedOnDefault = TRUE)
  
  ## test passes when results are rounded off to the second decimal place
  
  stopifnot(all.equal(round(result.1, 2), round(truth.1, 1)))
  
  
  ## results when we don't enter the rates manually are less accurate 
   
  result.2 <- upfront_to_spread(x = data.frame(date = as.Date("2014-04-22"),
                                               currency = "USD",
                                               coupon = 100,
                                               recovery = .4,
                                               tenor = 5,
                                               upfront = 18623.7),
                                isPriceClean = FALSE,
                                notional = 10000000,
                                payAccruedAtStart = FALSE,
                                payAccruedOnDefault = TRUE)
  
  stopifnot(all.equal(round(result.2), round(truth.1)))
  
  
  ## test cases to make sure results of the function don't change over time
  
  truth.2 <- upfront_to_spread(x = data.frame(date = as.Date("2014-01-14"),
                                              currency = "USD",
                                              coupon = 100,
                                              recovery = .4,
                                              tenor = 5,
                                              upfront = -10000000*3.48963/100),
                               isPriceClean = FALSE,
                               notional = 10000000,
                               payAccruedAtStart = FALSE,
                               payAccruedOnDefault = TRUE)
  
  truth.3 <- upfront_to_spread(x = data.frame(date = as.Date("2014-01-14"),
                                              currency = "USD",
                                              coupon = 100,
                                              recovery = .4,
                                              tenor = 5,
                                              upfront = -10000000*3.41/100),
                               isPriceClean = FALSE,
                               notional = 10000000,
                               payAccruedAtStart = TRUE,
                               payAccruedOnDefault = TRUE)
  
  result.3 <- upfront_to_spread(x = data.frame(date = as.Date("2014-01-14"),
                                               currency = "USD",
                                               coupon = 100,
                                               recovery = .4,
                                               tenor = 5,
                                               upfront = -10000000*3.48963/100),
                                isPriceClean = FALSE,
                                notional = 10000000,
                                payAccruedAtStart = FALSE,
                                payAccruedOnDefault = TRUE)
  
  stopifnot(all.equal(result.3, truth.2))
  
  
  result.4 <- upfront_to_spread(x = data.frame(date = as.Date("2014-01-14"),
                                               currency = "USD",
                                               coupon = 100,
                                               recovery = .4,
                                               tenor = 5,
                                               upfront = -10000000*3.41/100),
                                isPriceClean = FALSE,
                                notional = 10000000,
                                payAccruedAtStart = TRUE,
                                payAccruedOnDefault = TRUE)
  
  stopifnot(all.equal(result.4, truth.3))
  
})
