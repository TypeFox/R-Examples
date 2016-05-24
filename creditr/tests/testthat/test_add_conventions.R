context("Test illegal input to add_conventions")

test_that("error should occur if there is no currency.var ",{
  
  x <- data.frame(date = c(as.Date("2014-05-06"), as.Date("2014-05-07")))
  expect_error(add_conventions(x))
})

context("Test add_conventions")

test_that("test add conventions", {
  
  ## used independently
  
  x1 <- data.frame(date = c(as.Date("2014-05-06"), as.Date("2014-05-07")), currency = c("USD", "JPY"))
  result1 <- add_conventions(x1)
  
  ## joint usage with add_dates
  
  x2 <- data.frame(date = c(as.Date("2014-04-22"), as.Date("2014-04-22")),
                   currency = c("USD", "EUR"),
                   tenor = c(5, 5),
                   spread = c(120, 110),
                   coupon = c(100, 100),
                   recovery = c(0.4, 0.4),
                   notional = c(10000000, 10000000))
  
  x2 <- add_dates(x2)
  result2 <- add_conventions(x2)
})