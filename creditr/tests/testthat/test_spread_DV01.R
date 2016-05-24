context("Test spread_DV01")

## SpreadDV01.test.R using CDS data for Xerox Corp on April 22, 2014

test_that("test for spread_DV01", {

   truth.1 <- 4825.49

  x <- data.frame(date = as.Date("2014-04-22"),
                  currency = "USD",
                  tenor = 5,
                  spread = 105.8,
                  coupon = 100,
                  recovery = 0.4,
                  notional = 10000000,
                  stringsAsFactors = FALSE)
  
  result.1 <- spread_DV01(x)
  
  ## test case passes when results are rounded to the nearest whole number
  
  stopifnot(all.equal(round(result.1), round(truth.1)))
  
})
