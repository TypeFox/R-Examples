context("Test summary")

library(utils)

## use capture.output to test the summary() method

test_that("test summary() method", {
  
  cds <- CDS(date = as.Date("2014-05-07"), tenor = 5, spread = 50, coupon = 100)
  
  output <- capture.output(summary(cds))
  
  ## below expectations are hard coded, but are fine, because all we need to care about
  ## is that summary() still produces the same output after we redesign the class.
  
  expect_that(output[1], equals("Contract Type:                      SNAC   Date:                      2014-05-07"))
  
  expect_that(output[2], equals("Entity Name:                          NA   RED:                               NA"))
  
  expect_that(output[3], equals("Currency:                            USD   End Date:                  2019-06-20"))
  
  expect_that(output[4], equals("Spread:                               50   Coupon:                           100"))
  
  expect_that(output[5], equals("Upfront:                        -259,647   Spread DV01:                    5,022"))
  
  expect_that(output[6], equals("IR DV01:                           64.44   Rec Risk (1 pct):               87.88"))
  }
)