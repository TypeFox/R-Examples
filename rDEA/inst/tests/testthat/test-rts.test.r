
data("hospitals", package="rDEA")

## choosing inputs and outputs for analysis
firms = 1:40
X = hospitals[firms, c('labor', 'capital')]
Y = hospitals[firms, c('inpatients', 'outpatients')]
W = hospitals[firms, c('labor_price', 'capital_price')]

B = 20

context("Input rts.test")

test_that("input rts.test with H0='constant'", {
  r = rts.test(X=X, Y=Y, model="input", H0="constant", B = B)
  expect_equal( length(r$theta_H0_hat), length(firms) )
  expect_equal( length(r$theta_vrs_hat), length(firms) )
  expect_equal( length(r$w_hat_boot),   B )
  expect_equal( length(r$w48_hat_boot), B )
  
  expect_true( all(r$theta_H0_hat <= 1) )
  expect_true( all(r$theta_vrs_hat <= 1) )
})

test_that("input rts.test with H0='non-increasing'", {
  r = rts.test(X=X, Y=Y, model="input", H0="non-increasing", B = B)
  expect_equal( length(r$theta_H0_hat), length(firms) )
  expect_equal( length(r$theta_vrs_hat), length(firms) )
  expect_equal( length(r$w_hat_boot),   B )
  expect_equal( length(r$w48_hat_boot), B )
  
  expect_true( all(r$theta_H0_hat <= 1) )
  expect_true( all(r$theta_vrs_hat <= 1) )
})

context("Output rts.test")

test_that("output rts.test with H0='constant'", {
  r = rts.test(X=X, Y=Y, model="output", H0="constant", B = B)
  expect_equal( length(r$theta_H0_hat), length(firms) )
  expect_equal( length(r$theta_vrs_hat), length(firms) )
  expect_equal( length(r$w_hat_boot),   B )
  expect_equal( length(r$w48_hat_boot), B )
  
  expect_true( all(r$theta_H0_hat <= 1) )
  expect_true( all(r$theta_vrs_hat <= 1) )
})

test_that("output rts.test with H0='non-increasing'", {
  r = rts.test(X=X, Y=Y, model="output", H0="non-increasing", B = B)
  expect_equal( length(r$theta_H0_hat), length(firms) )
  expect_equal( length(r$theta_vrs_hat), length(firms) )
  expect_equal( length(r$w_hat_boot),   B )
  expect_equal( length(r$w48_hat_boot), B )
  
  expect_true( all(r$theta_H0_hat <= 1) )
  expect_true( all(r$theta_vrs_hat <= 1) )
})
