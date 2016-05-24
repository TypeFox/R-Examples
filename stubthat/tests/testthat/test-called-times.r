library(testthat)

simpf <- function(a = 1, b, d, ...) return(5)
not_expected_error <- 'Function is called with arguments different from expected!'

test_that('calledTimes: returns the number of times the function has been called', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$returns(10)
  
  temp <- stub_of_simpf$f(1, 2, 3, 4)
  temp <- stub_of_simpf$f(1, 3, 3, b = 4)
  
  expect_equal(stub_of_simpf$calledTimes(), 2)
  
  temp <- stub_of_simpf$f(1, 2, 3, 4)
  
  expect_equal(stub_of_simpf$calledTimes(), 3)
})