context("system utils")

test_that("sysutils",{
  expect_is(n<-ncpus(),'integer')
  expect_true(n>=1L)
})
