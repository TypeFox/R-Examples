library(easyVerification)
context('Unwrapping')

xxx <- array(1, c(1,3,5))
xx <- array(1, c(3,5))
x <- rep(1, 15)

test_that("input arguments are valid", {
  expect_error(veriUnwrap(xxx, 'EnsMe'))
  expect_error(veriUnwrap(x, 'EnsMe'))
  expect_error(veriUnwrap(xx, 'fdkl'))
  expect_error(veriUnwrap(xx, EnsMe))
})


test_that("returned value is correct", {
  expect_is(veriUnwrap(xx, 'EnsMe'), 'numeric')
  expect_is(veriUnwrap(xx, 'FairCrps'), mode(FairCrps(xx[,-ncol(xx)], xx[,ncol(xx)])))
  expect_equal(veriUnwrap(xx, 'EnsMe'), 0)
})

