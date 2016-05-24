context("rToJags")

test_that("rToJags: exponents",
{
  expect_equal(rToJags(y ~ x^2),
               "y ~ pow(x,2)")
})

test_that("rToJags: logit (from VGAM package)",
{
  expect_equal(rToJags(y ~ logit(x, inverse=TRUE)),
               "y ~ ilogit(x)")
})