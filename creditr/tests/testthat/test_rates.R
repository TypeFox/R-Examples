context("Test values in rates data frame")

data(rates)

test_that("test that rates data frame has correct variable names and types", {
  
  expect_equal(names(rates), c("date", "currency", "expiry", "rate"))
  expect_equal(unname(sapply(rates, class)), c("Date", "character", "character", "numeric"))  

})

test_that("test that rates data frame variables currency and expiry have only allowed values", {
  
  expect_true(all(rates$currency %in% c("USD", "JPY", "EUR")))
  
  expect_true(all(rates$expiry %in% c("1M", "2M", "3M", "6M", "9M",
                                      "1Y", "2Y", "3Y", "4Y", "5Y",  
                                      "6Y", "7Y", "8Y", "9Y",
                                      "10Y", "12Y", "15Y", "20Y", "25Y", "30Y")))
  
  
})

## independence day test case

test_that("test that holidays are covered in rates.RData", {
  
  ## if we trade on 2012-07-05, then we should use the previous business day's
  ## interest rate, which is 2012-07-03. But notice here, that the interest
  ## rate of 2012-07-03 is not stored in rates.RData as the row with date
  ## 2012-07-03; instead, it's the row with 2012-07-04 that contains the 
  ## true interest rate of 2012-07-03, because we have to adjust one
  ## business day after 2012-07-03 in rates.RData
  
  rate.holi.1 <- rates[rates$currency == "USD" & rates$date == as.Date("2012-07-05") & 
                       rates$expiry == "1M", ]$rate
  rate.holi.2 <- rates[rates$currency == "USD" & rates$date == as.Date("2012-07-04") & 
                         rates$expiry == "1M", ]$rate
  expect_equal(rate.holi.1, rate.holi.2)
})

## a random weekend test case. Done for each currency separately, 
## on separate week-ends.

test_that("test that weekends are covered correctly", {
  
  ## a random weekend for USD
  
  ## we want to see that for sunday, 2014-08-10, its interest rate should be same with
  ## saturday, 2014-08-09, because their previous business days are both
  ## friday, 2014-08-08.
  
  rate.weekend.1 <- rates[rates$currency == "USD" & rates$date == as.Date("2012-08-12") & 
                            rates$expiry == "1M", ]$rate
  rate.weekend.2 <- rates[rates$currency == "USD" & rates$date == as.Date("2012-08-11") & 
                            rates$expiry == "1M", ]$rate
  expect_equal(rate.weekend.1, rate.weekend.2)
  ## a random weekend for JPY
  
  rate.weekend.3 <- rates[rates$currency == "JPY" & rates$date == as.Date("2012-08-05") & 
                            rates$expiry == "1M", ]$rate
  rate.weekend.4 <- rates[rates$currency == "JPY" & rates$date == as.Date("2012-08-04") & 
                            rates$expiry == "1M", ]$rate
  expect_equal(rate.weekend.3, rate.weekend.4)
  ## a random weekend for EUR
  
  rate.weekend.5 <- rates[rates$currency == "EUR" & rates$date == as.Date("2012-08-25") & 
                            rates$expiry == "1M", ]$rate
  rate.weekend.6 <- rates[rates$currency == "EUR" & rates$date == as.Date("2012-08-26") & 
                            rates$expiry == "1M", ]$rate
  expect_equal(rate.weekend.5, rate.weekend.6)
  
})


test_that("Test to show there are no missing dates", {
  
  expect_equal(length(unique(rates$date)), 
               as.numeric(max(rates$date) - min(rates$date) + 1))
  
})


test_that("test that rates don't move `too much' day-over-day", {

  ## Problem with this test is that what is `too much' movement depends on the
  ## expiry. 
    
  sample.df <- rates[rates$currency == "EUR" & rates$expiry == "20Y",]
    
  sample.date <- sample(sample.df$date, size = 1)
  sample.date.next <- sample.date + 1
    
  rates.1 <- rates[rates$date == sample.date & rates$currency == "EUR" & 
                       rates$expiry == "20Y", ]$rate
    
  rates.2 <- rates[rates$date == sample.date.next & rates$currency == "EUR" & 
                       rates$expiry == "20Y", ]$rate
    
  expect_true(abs(rates.2-rates.1)/rates.1 < 0.1)
})

