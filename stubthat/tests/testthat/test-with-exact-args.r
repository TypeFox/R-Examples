library(testthat)

simpf <- function(a = 1, b, d, ...) return(5)
not_expected_error <- 'Function is called with arguments different from expected!'

test_that('Returns the specified value if called with the exact arguments specified', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withExactArgs(a = 1, b = 2, d = 3, c = 4)$returns(10)
  stub_func <- stub_of_simpf$f
  
  expect_equal(stub_func(1, 2, 3, c = 4), 10)
  expect_null(stub_func(2, 2, 3))
  expect_equal(stub_func(1, 2, 3, c = 4), 10)
  expect_equal(stub_func(1, 2, c = 4, 3), 10)
})

test_that('Throws error with specified message if called with the exact arguments specified', {
  stub_of_simpf <- stub(simpf)
  err_msg <- 'error is good'
  stub_of_simpf$withExactArgs(a = 1, b = 2, d = 3, c = 4)$throws(err_msg)
  stub_func <- stub_of_simpf$f
  
  expect_error(stub_func(1, 2, 3, c = 4), err_msg)
  expect_null(stub_func(2, 2, 3))
  expect_error(stub_func(1, 2, 3, c = 4), err_msg)
  expect_error(stub_func(1, c = 4, 2, 3), err_msg)
})

test_that('It does the right thing even when there are multiple expectations - withExactArgs.return/throw and default return', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withExactArgs(a = 1, b = 2, d = 3, c = 4)$returns(10)
  stub_of_simpf$withExactArgs(a = 2, b = 2, d = 3, c = 4)$throws('err')
  stub_of_simpf$returns('a')
  stub_func <- stub_of_simpf$f

  expect_equal(stub_func(1, 2, 3, c = 4), 10)
  expect_equal(stub_func(3, 2, 3, c = 4), 'a')
  expect_error(stub_func(2, 2, 3, c = 4), 'err')
  
  expect_equal(stub_func(3, 2, 3, c = 4), 'a')
  expect_error(stub_func(2, c = 4, 2, 3), 'err')
  expect_equal(stub_func(c = 4, 1, 2, 3), 10)
})

test_that('It does the right thing even when there are multiple expectations - withExactArgs.return/throw and default throw', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withExactArgs(a = 1, b = 2, d = 3, c = 4)$returns(10)
  stub_of_simpf$withExactArgs(a = 2, b = 2, d = 3, c = 4)$throws('err')
  stub_of_simpf$throws('pqrs')
  stub_func <- stub_of_simpf$f

  expect_equal(stub_func(1, 2, 3, c = 4), 10)
  expect_error(stub_func(3, 2, 3, c = 4), 'pqrs')
  expect_error(stub_func(2, 2, 3, c = 4), 'err')
  
  expect_error(stub_func(c = 4, 3, 2, 3), 'pqrs')
  expect_equal(stub_func(1, 2, 3, c = 4), 10)
  expect_error(stub_func(2, 2, 3, c = 4), 'err')
})
