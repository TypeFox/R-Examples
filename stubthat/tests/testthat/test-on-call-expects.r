library(testthat)

simpf <- function(a = 1, b, d, ...) return(5)
not_expected_error <- 'Function is called with arguments different from expected!'

test_that('expects: It returns the specified value when expected arguments are part of the function call on the nth time running of the function', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$onCall(3)$expects(g = 'a')$returns(10)
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_equal(stub_func(1, 2, 3, g = 'a'), 10)
  expect_null(stub_func(1, 2, 3, g = 'a'))
})

test_that('expects: It throws error when expected arguments are not part of the function call on the nth call', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$onCall(3)$expects(g = 'a')$returns(10)
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_error(stub_func(1, 2, 3, g = 'b'), not_expected_error)
  expect_null(stub_func(1, 2, 3, g = 'a'))
})

test_that('expects: It throws error with the specified message when expected arguments are part of the function call on the nth time running of the function', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$onCall(3)$expects(g = 'a')$throws('error is nice')
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_error(stub_func(1, 2, 3, g = 'a'), 'error is nice')
  expect_null(stub_func(1, 2, 3, g = 'a'))
})

test_that('expects: It throws error when expected arguments are not part of the function call on the nth call', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$onCall(3)$expects(g = 'a')$throws('error is nice')
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_null(stub_func(1, 2, 3, g = 'a'))
  expect_error(stub_func(1, 2, 3, g = 'b'), not_expected_error)
  expect_null(stub_func(1, 2, 3, g = 'a'))
})
