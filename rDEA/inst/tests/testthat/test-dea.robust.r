
library("maxLik")

data("hospitals", package="rDEA")

## choosing inputs and outputs for analysis
firms = 1:50
X = hospitals[firms, c('labor', 'capital')]
Y = hospitals[firms, c('inpatients', 'outpatients')]
W = hospitals[firms, c('labor_price', 'capital_price')]

context("Input DEA robust")

test_that("input DEA robust with variable RTS", {
  dr = dea.robust(X=X, Y=Y, model="input", RTS="variable", B=20)
  expect_equal( length(dr$bias),          length(firms) )
  expect_equal( length(dr$theta_hat),     length(firms) )
  expect_equal( length(dr$theta_hat_hat), length(firms) )
  expect_equal( length(dr$theta_ci_low),  length(firms) )
  expect_equal( length(dr$theta_ci_high), length(firms) )
  expect_equal( nrow(dr$theta_hat_star),  length(firms) )
  expect_equal( ncol(dr$theta_hat_star),  20 )
  
  expect_true( all(dr$bias >= 0) )
  expect_true( all(dr$ci_low <= dr$ci_high) )
})

test_that("input DEA robust with variable RTS and bw_mult=1.5", {
  dr = dea.robust(X=X, Y=Y, model="input", RTS="variable", B=20, bw_mult = 1.5)
  expect_equal( length(dr$bias),          length(firms) )
  expect_equal( length(dr$theta_hat),     length(firms) )
  expect_equal( length(dr$theta_hat_hat), length(firms) )
  expect_equal( length(dr$theta_ci_low),  length(firms) )
  expect_equal( length(dr$theta_ci_high), length(firms) )
  expect_equal( nrow(dr$theta_hat_star),  length(firms) )
  expect_equal( ncol(dr$theta_hat_star),  20 )
  
  expect_true( all(dr$bias >= 0) )
  expect_true( all(dr$ci_low <= dr$ci_high) )
})

context("Output DEA robust")

test_that("output DEA robust with variable RTS", {
  dr = dea.robust(X=X, Y=Y, model="output", RTS="variable", B=20)
  expect_equal( length(dr$bias),          length(firms) )
  expect_equal( length(dr$theta_hat),     length(firms) )
  expect_equal( length(dr$theta_hat_hat), length(firms) )
  expect_equal( length(dr$theta_ci_low),  length(firms) )
  expect_equal( length(dr$theta_ci_high), length(firms) )
  expect_equal( nrow(dr$theta_hat_star),  length(firms) )
  expect_equal( ncol(dr$theta_hat_star),  20 )
  
  expect_true( all(dr$bias >= 0) )
  expect_true( all(dr$ci_low <= dr$ci_high) )
})

context("Costmin(Fare) DEA robust")

test_that("Costmin(Fare) DEA robust with variable RTS", {
  dr = dea.robust(X=X, Y=Y, W=W, model="costmin", RTS="variable", B=20)
  expect_equal( length(dr$bias),          length(firms) )
  expect_equal( length(dr$gamma_hat),     length(firms) )
  expect_equal( length(dr$gamma_hat_hat), length(firms) )
  expect_equal( length(dr$gamma_ci_low),  length(firms) )
  expect_equal( length(dr$gamma_ci_high), length(firms) )
  expect_equal( nrow(dr$gamma_hat_star),  length(firms) )
  expect_equal( ncol(dr$gamma_hat_star),  20 )
  
  expect_true( all(dr$bias >= 0) )
  expect_true( all(dr$gamma_ci_low <= dr$gamma_ci_high) )
})

