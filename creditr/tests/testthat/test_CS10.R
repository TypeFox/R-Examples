context("Test CS10")

test_that("test case for CS10", {
  x <- data.frame(date = c(as.Date("2014-04-22"), as.Date("2014-04-22")),
                  currency = c("USD", "EUR"),
                  tenor = c(5, 5),
                  spread = c(105.8, 99),
                  coupon = c(100, 100),
                  recovery = c(0.4, 0.4),
                  notional = c(10000000, 10000000),
                  stringsAsFactors = FALSE)
  
  result <- CS10(x)
  ## we don't have any thing to test this against at the moment
})