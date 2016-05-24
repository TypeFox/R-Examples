context("Testing information_arma")


test_that("bogus arguments throw error",{
  expect_error(information_arma(NA, Inf))
})

test_that("output of information_arma is correct",{
  mat <- information_arma(c(0.9, 0.45))
  expect_is(mat, "matrix")
  expect_equal(mat, matrix(-c(0.7474095,1.2230338,1.2230338,0.7474095), 2, 2), tolerance = 1e-5)

  mat <- information_arma(theta = c(0.9, 0.45))
  expect_is(mat, "matrix")
  expect_equal(mat, matrix(c(2.039740,-1.266045,-1.266045,2.039740), 2, 2), tolerance = 1e-5)

  mat <- information_arma(0.9,0.45)
  expect_is(mat, "matrix")
  expect_equal(mat, matrix(c(5.2631579,-0.7117438,-0.7117438,1.2539185), 2, 2), tolerance = 1e-5)
})
