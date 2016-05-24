context("Test show")

library(utils)

## use capture.output to test the show() method

test_that("test show() method", {
  
  cds <- CDS(date = as.Date("2014-05-07"), tenor = 5, spread = 50, coupon = 100)
  
  output <- capture.output(show(cds))
  
  expect_that(output[1], equals("CDS Contract "))
  
  expect_that(output[2], equals("Contract Type:                      SNAC   Currency:                         USD"))
  
  expect_that(output[3], equals("Entity Name:                          NA   RED:                               NA"))
  
  expect_that(output[4], equals("date:                         2014-05-07"))
  
  expect_that(output[5], equals(""))
  
  expect_that(output[6], equals("Calculation "))
  
  expect_that(output[7], equals("price:                            102.46   Spread:                            50"))
  
  expect_that(output[8], equals("Principal:                      -246,036   Spread DV01:                    5,022"))
  
  ## Actually, hard coding doesn't matter here since we are dealing with this
  ## specific case.
  
  ## All we have to make sure is that after we redesign the class these prints
  ## still look the same
  
  expect_that(output[9], equals("Accrual:                         -13,611   IR DV01:                        64.44"))
  
  expect_that(output[10], equals("Upfront:                        -259,647   Rec Risk (1 pct):               87.88"))
  
  expect_that(output[11], equals("Default Prob:                     0.0424"))
  
  expect_that(output[12], equals(""))}
)
