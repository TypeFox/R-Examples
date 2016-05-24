library(testthat)

simpf <- function(a = 1, b, d, ...) return(5)
not_expected_error <- 'Function is called with arguments different from expected!'

test_that('strictlyExpects: It returns the specified value when called with the exact set of expected arguments on the nth time running of the function', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$onCall(3)$strictlyExpects(a = 1, b = 2, d = 3, c = 4)$returns(10)
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, c = 4))
  expect_null(stub_func(1, 2, 3, c = 4))
  expect_equal(stub_func(1, 2, 3, c = 4), 10)
  expect_null(stub_func(1, 2, 3, c = 4))
})

test_that('strictlyExpects: It throws error when not called with the exact set of expected arguments on the nth call', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$onCall(3)$strictlyExpects(a = 1, b = 2, d = 3, c = 4)$returns(10)
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, c = 4))
  expect_null(stub_func(1, 2, 3, c = 4))
  expect_error(stub_func(1, 2, 3, c = 5), not_expected_error)
  expect_null(stub_func(1, 2, 3, c = 4))
})

test_that('strictlyExpects: It throws error with the specified message when called with the exact set of expected arguments on the nth time running of the function', {
  stub_of_simpf <- stub(simpf)
  err_msg <- 'error is good'
  stub_of_simpf$onCall(2)$strictlyExpects(a = 1, b = 2, d = 3, c = 4)$throws(err_msg)
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, c = 4))
  expect_error(stub_func(1, 2, 3, c = 4), err_msg)
  expect_null(stub_func(1, 2, 3, c = 4))
})

test_that('strictlyExpects: It throws error when not called with the exact set of expected arguments on the nth call', {
  stub_of_simpf <- stub(simpf)
  err_msg <- 'error is good'
  stub_of_simpf$onCall(2)$strictlyExpects(a = 1, b = 2, d = 3, c = 4)$throws(err_msg)
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, c = 4))
  expect_error(stub_func(1, 2, 3, c = 5), not_expected_error)
  expect_null(stub_func(1, 2, 3, c = 4))
})

