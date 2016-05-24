context("Test IR_DV01")

test_that("test for IR_DV01", {
  
  ## comparing IR.DV01 calculated by our package for Xerox Corp and Electrolux
  ## AB on April 22, 2014 with the results on Bloomberg
  
  x <- data.frame(date = c(as.Date("2014-04-22"), as.Date("2014-04-22")),
                  currency = c("USD", "EUR"),
                  tenor = c(5, 5),
                  spread = c(105.8, 99),
                  coupon = c(100, 100),
                  recovery = c(0.4, 0.4),
                  notional = c(10000000, 10000000),
                  stringsAsFactors = FALSE)
  
  result <- IR_DV01(x)
  
  truth <- c(-7.36, 1.29)
  
  stopifnot(all.equal(round(result), round(truth)))
})
