context("set_pars")

test_that("check length of output",{
  expect_equal(length(.set_pars("Bailey2001")), 12)
  expect_equal(length(.set_pars("Bailey2002")), 12)
  expect_equal(length(.set_pars("Bailey2004")), 12)
  expect_equal(length(.set_pars("Pagonis2007")), 12)
  expect_equal(length(.set_pars("Pagonis2008")), 12)
})
