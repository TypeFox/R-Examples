context("Test check_inputs") 

test_that("test for check.input", {
  
  ## define a valid data frame of inputs. The test below just modifies it.
  
  x0 <- data.frame(dates = c(as.Date("2014-04-22"), as.Date("2014-04-22")),
                   currency = c("USD", "EUR"),
                   maturity = c(as.Date("2019-04-22"), as.Date("2019-04-22")),
                   tenor = c(5, 5),
                   spread = c(120, 110),
                   coupon = c(100, 100),
                   recovery = c(0.4, 0.4),
                   notional = c(10000000, 10000000))
  
  ## if inputs are missing, it should return an error
  
  x1 <- data.frame(dates = c(as.Date("2014-04-22"), as.Date("2014-04-22")))
  expect_error(check_inputs(x1))
  
  ## if date is not a Date class
  x2 <- x0
  x2$dates <- c(20140102, 20140309)
  expect_error(check_inputs(x2))
  
  ## if tenor is not a numeric class
  
  
  x3 <- x0
  x3$tenor <- c("5Y", "5Y")
  expect_error(check_inputs(x3))
  
  
  ## if maturity is not a Date class
  
  x4 <- x0
  x4$maturity <- c(20190101, 20190102)
  expect_error(check_inputs(x4))
  
  ## if spread is not a numeric class
  
  x5 <- x0
  x5$spread <- c("500", "500")
  expect_error(check_inputs(x5))
  
  ## if coupon is not a numeric class
  
  x6 <- x0
  x6$coupon <- c("100", "100")
  expect_error(check_inputs(x6))
  
  ## if notional is not a numeric class
  
  x7 <- x0
  x7$notional <- c("10000000", "10000000")
  expect_error(check_inputs(x7))
})
