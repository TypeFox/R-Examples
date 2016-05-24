library(stubthat)
library(testthat)

simpf <- function(a = 1, d, ...) return(5)
not_expected_error <- 'Function is called with arguments different from expected!'

test_that('Testing with function given as input', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$expects(d = paste0)
  stub_func <- stub_of_simpf$f
  expect_null(stub_func(d = paste0))
  
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$strictlyExpects(a = 1, d = paste0, e = 'a')
  stub_of_simpf$onCall(2)$expects(d = paste)
  stub_of_simpf$onCall(3)$strictlyExpects(d = paste)
  stub_func <- stub_of_simpf$f
  expect_null(stub_func(d = paste0, e = 'a'))
  expect_null(stub_func(paste, a = 10))
  expect_error(stub_func(paste, a = 10))
})
