context("Testing acv_arma")

test_that("bogus arguments throw error",{
  expect_error(acv_arma(Inf,1,10))
})

test_that("output of acv_arma is correct",{
  n <- 0:9
  answer <- 2^(-n) * (32/3 + 8 * n) /(32/3)
  acv <- acv_arma(c(1.0, -0.25), 1.0, 10)
  expect_equal(acv/acv[1], answer)
})
