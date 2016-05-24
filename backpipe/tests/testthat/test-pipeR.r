library(pipeR)

context("backpipe with pipeR")

test_that("backpipe works with pipeR as intended", {
  
  expect_equal( mean %<<% 1:3, 2 )
  expect_equal( mean %<<% range( na.rm = TRUE ) %<<% c(1:3, NA), 2 )
  
})
