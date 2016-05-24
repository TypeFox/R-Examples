
library("maxLik")
data("hospitals", package="rDEA")

## choosing inputs and outputs for analysis
firms = 1:40
X = hospitals[firms, c('labor', 'capital')]
Y = hospitals[firms, c('inpatients', 'outpatients')]
W = hospitals[firms, c('labor_price', 'capital_price')]
Z = hospitals[firms, c('z1'), drop = FALSE]

L1 = 20
L2 = 20
Nbeta = ncol(Z) + 1

context("Input DEA Environment robust")

test_that("input DEA robust with variable RTS", {
  der = dea.env.robust(X=X, Y=Y, Z=Z, model="input", RTS="variable", L1=L1, L2=L2)
  expect_equal( length(der$bias),          length(firms) )
  expect_equal( length(der$delta_hat),     length(firms) )
  expect_equal( length(der$delta_hat_hat), length(firms) )
  expect_equal( length(der$delta_ci_low),  length(firms) )
  expect_equal( length(der$delta_ci_high), length(firms) )
  expect_equal( length(der$delta_ci_kneip_low),  length(firms) )
  expect_equal( length(der$delta_ci_kneip_high), length(firms) )
  expect_equal( length(der$beta_hat),          Nbeta )
  expect_equal( length(der$beta_hat_hat),      Nbeta )
  expect_equal( length(der$beta_hat_hat_star), Nbeta )
  
  expect_true( all(der$bias <= 0) )
  expect_true( all(der$delta_hat >= 1) )
})

context("Output DEA Environment robust")

test_that("output DEA robust with variable RTS", {
  der = dea.env.robust(X=X, Y=Y, Z=Z, model="output", RTS="variable", L1=L1, L2=L2)
  expect_equal( length(der$bias),          length(firms) )
  expect_equal( length(der$delta_hat),     length(firms) )
  expect_equal( length(der$delta_hat_hat), length(firms) )
  expect_equal( length(der$beta_hat),          Nbeta )
  expect_equal( length(der$beta_hat_hat),      Nbeta )
  expect_equal( length(der$beta_hat_hat_star), Nbeta )
  
  expect_true( all(der$bias <= 0) )
  expect_true( all(der$delta_hat >= 1) )
})

context("Costmin(Fare) DEA Environment robust")

test_that("costmin(fare) DEA robust with variable RTS", {
  der = dea.env.robust(X=X, Y=Y, W=W, Z=Z, model="costmin", RTS="variable", L1=L1, L2=L2)
  expect_equal( length(der$bias),          length(firms) )
  expect_equal( length(der$delta_hat),     length(firms) )
  expect_equal( length(der$delta_hat_hat), length(firms) )
  expect_equal( length(der$beta_hat),          Nbeta )
  expect_equal( length(der$beta_hat_hat),      Nbeta )
  expect_equal( length(der$beta_hat_hat_star), Nbeta )
  
  expect_true( all(der$bias <= 0) )
  expect_true( all(der$delta_hat >= 1) )
})

