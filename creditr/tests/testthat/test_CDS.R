context("Test CDS")

## CDS.R test case for Caesars Entertainment Operating Co Inc using data obtainined from Bloomberg.
## The .png files containing these results are in tests/sources

## rel.diff is a function which determines if the relative difference
## between the first argument (truth) and the second (calculated) is within
## a designated acceptable range (threshold).

rel.diff <- function(truth,
                     calculated,
                     threshold = 0.01){
  return (abs(calculated - truth) / truth < threshold)
}

## we use rates for the relevant date and currency, extracted from the rates data frame, stored in the data 
## directory. 

test_that("test for the CDS", {
  
  result.1 <- CDS(date = as.Date("2014-04-15"),
                  currency = "USD",
                  maturity = as.Date("2019-06-20"),                    
                  spread = 12354.529,
                  coupon = 500,
                  recovery = 0.4,
                  notional = 10000000)
  
  ## The following result matches perfectly with Bloomberg
  
  expect_equal(42.55, round(result.1@price, 2))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
 
  expect_true(rel.diff(5707438, result.1@upfront))
  
  expect_true(rel.diff(-271.18, result.1@IR.DV01))
  
  expect_true(rel.diff(5744938, result.1@principal))
  
  expect_true(rel.diff(-95430.32, result.1@rec.risk.01))
  
  expect_true(rel.diff(21.15, result.1@spread.DV01))
  
  
  ## CDS.R test case for Chorus Ltd. (Australian company)
  
  result.2 <- CDS(date = as.Date("2014-04-15"),
                  currency = "USD",                    
                  maturity = as.Date("2019-06-20"),                    
                  spread = 243.28,
                  coupon = 100,
                  recovery = 0.40,
                  notional = 10000000)
  
  ## The following results match perfectly with Bloomberg
  
  expect_equal(650580, round(result.2@upfront))
  
  expect_equal(658080, round(result.2@principal))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
  
  expect_true(rel.diff(-169.33, result.2@IR.DV01))
  
  expect_true(rel.diff(93.42, result.2@price))
  
  expect_true(rel.diff(-1106.34, result.2@rec.risk.01))
  
  expect_true(rel.diff(4317.54, result.2@spread.DV01))

  
  
  ## CDS.R test case for Electrolux AB corporation
  
  result.3 <- CDS(date = as.Date("2014-04-22"),
                  tenor = 5,
                  spread = 99,
                  contract ="STEC",
                  currency="EUR",
                  coupon = 100,
                  recovery = 0.4,
                  notional = 10000000)
  
  ## The following results match perfectly with Bloomberg
  
  expect_equal(-14368, round(result.3@upfront))
  
  expect_equal(-4924, round(result.3@principal))
  
  expect_equal(100.05, round(result.3@price, 2))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
  
  expect_true(rel.diff(1.29, result.3@IR.DV01))
  
  expect_true(rel.diff(3.46, result.3@rec.risk.01, threshold = 0.1))

  expect_true(rel.diff(4923.93, result.3@spread.DV01))

  
  ## CDS.R test case for Norske Skogindustrier ASA (European company)
  
  result.4 <- CDS(date = as.Date("2014-04-15"),
                  tenor = 5,
                  spread = 2785.8889,
                  contract ="STEC",
                  currency="EUR",
                  coupon = 500,
                  recovery = 0.4,
                  notional = 10000000)
  
  ## The following results match perfectly with Bloomberg
  
  expect_equal(55.5, round(result.4@price, 1))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
  
  expect_true(rel.diff(4412500, result.4@upfront))
  
  expect_true(rel.diff(-727.47, result.4@IR.DV01))

  expect_true(rel.diff(4450000, result.4@principal))
  
  expect_true(rel.diff(-56413.77, result.4@rec.risk.01))
  
  expect_true(rel.diff(731.48, result.4@spread.DV01))

  
  ## CDS.R test case for RadioShack Corp
  
  result.5 <- CDS(date = as.Date("2014-04-15"),
                  currency = "USD",                     
                  maturity = as.Date("2019-06-20"),
                  spread = 9106.8084,
                  coupon = 500,
                  recovery = 0.4,
                  notional = 10000000)
  
  ## The following results match perfectly with Bloomberg
  
  expect_equal(43.5, round(result.5@price, 1))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
  
  expect_true(rel.diff(5612324, result.5@upfront))
  
  expect_true(rel.diff(-361.62, result.5@IR.DV01))

  expect_true(rel.diff(5649824, result.5@principal))
  
  expect_true(rel.diff(-93430.52, result.5@rec.risk.01))
  
  expect_true(rel.diff(40.86, result.5@spread.DV01))
  
  
  ## CDS.R test case for Tokyo Electric Power Co. Inc.
  
  result.6 <- CDS(date = as.Date("2014-04-15"),
                  tenor = 5,
                  contract ="STEC",
                  spread = 250,
                  currency = "JPY",
                  coupon = 100,
                  recovery = 0.35,
                  notional = 10000000)
  
  ## The following results match perfectly with Bloomberg
  
  expect_equal(92.91, round(result.6@price, 2))
  
  expect_equal(701502, round(result.6@upfront))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
  
  expect_true(rel.diff(-184.69, result.6@IR.DV01))

  expect_true(rel.diff(709002, result.6@principal))
  
  expect_true(rel.diff(-1061.74, result.6@rec.risk.01))
  
  expect_true(rel.diff(4448.92, result.6@spread.DV01))
  
  
  ## CDS.R test case for Toys R Us Inc
  
  result.7 <- CDS(date = as.Date("2014-04-15"),
                  tenor = 5,
                  contract="SNAC",
                  spread = 1737.7289,
                  currency = "USD",
                  coupon = 500,
                  recovery = 0.40,
                  notional = 10000000)
  
  ## The following results match perfectly with Bloomberg
  
  expect_equal(3237500, round(result.7@upfront))
  
  expect_equal(3275000, round(result.7@principal))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
  
  expect_true(rel.diff(-648.12, result.7@IR.DV01))
  
  expect_true(rel.diff(67.25, result.7@price))

  expect_true(rel.diff(-30848.67, result.7@rec.risk.01))

  expect_true(rel.diff(1580.31, result.7@spread.DV01))
  
  ## CDS.R test case for Xerox corporation
  
  result.8 <- CDS(date = as.Date("2014-04-22"),
                  tenor = 5,
                  spread = 105.8,
                  coupon = 100,
                  recovery = 0.4,
                  notional = 10000000)
  
  ## The following results match perfectly with Bloomberg
  
  expect_equal(18624, round(result.8@upfront))
  
  expect_equal(28068, round(result.8@principal))
  
  ## The following results do not match exactly, so we will check
  ## to see if our calcualted results are within an acceptable range
  ## of the true values from Bloomberg
  
  expect_true(rel.diff(-7.36, result.8@IR.DV01))
  
  expect_true(rel.diff(99.71931785, result.8@price))
  
  expect_true(rel.diff(-20.85, result.8@rec.risk.01))

  expect_true(rel.diff(4825.49, result.8@spread.DV01))}
)