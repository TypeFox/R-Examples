library(testthat)

simpf <- function(a = 1, b, d, ...) return(5)
not_expected_error <- 'Function is called with arguments different from expected!'

test_that('returns: Always returns specified value', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$returns(10)
  stub_func <- stub_of_simpf$f
  expect_equal(stub_func(1, 2, 3, 4), 10)
  expect_equal(stub_func(1, 3, 3, b = 4), 10)
  expect_equal(stub_func(1, a = 5, 3, 4), 10)
})

test_that('throws: Always throws error with specified msg', {
  stub_of_simpf <- stub(simpf)
  throw_msg <- 'err msg xyz'
  stub_of_simpf$throws(throw_msg)
  stub_func <- stub_of_simpf$f
  expect_error(stub_func(1, 2, 3, 4), throw_msg)
  expect_error(stub_func(1, 3, 3, b = 4), throw_msg)
  expect_error(stub_func(1, a = 5, 3, 4), throw_msg)
})

test_that('strictlyExpects: Always checks the function call with expected arguments (exact set)', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$strictlyExpects(a = 1, b = 2, d = 3, c = 4)
  stub_func <- stub_of_simpf$f

  expect_null(stub_func(1, 2, 3, c = 4))
  expect_null(stub_func(c = 4, 2, a = 1, 3))

  expect_error(stub_func(2, 3, 3, c = 4), not_expected_error)
  expect_error(stub_func(2, 3, 3), not_expected_error)
  expect_error(stub_func(c = 4, a = 3, 1, 2), not_expected_error)
  expect_error(stub_func(a = 3, 1, 2), not_expected_error)
})

test_that('strictlyExpects & returns: Always checks the function call with expected arguments (exact set) and returns the specified value', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$strictlyExpects(a = 1, b = 2, d = 3, c = 4)
  stub_of_simpf$returns('a')
  stub_func <- stub_of_simpf$f

  expect_equal(stub_func(1, 2, 3, c = 4), 'a')
  expect_equal(stub_func(c = 4, 2, a = 1, 3), 'a')

  expect_error(stub_func(2, 3, 3, c = 4), not_expected_error)
  expect_error(stub_func(2, 3, 3), not_expected_error)
  expect_error(stub_func(c = 4, a = 3, 1, 2), not_expected_error)
})

test_that('strictlyExpects & throws: Always checks the function call with expected arguments (exact set) and throws error with specified msg', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$strictlyExpects(a = 1, b = 2, d = 3, c = 4)
  stub_of_simpf$throws('err msg')
  stub_func <- stub_of_simpf$f

  expect_error(stub_func(1, 2, 3, c = 4), 'err msg')
  expect_error(stub_func(c = 4, 2, a = 1, 3), 'err msg')

  expect_error(stub_func(2, 3, 3, c = 4), not_expected_error)
  expect_error(stub_func(2, 3, 3), not_expected_error)
  expect_error(stub_func(c = 4, a = 3, 1, 2), not_expected_error)
})

test_that('expects: Always checks if the expected arguments are part of the function call', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$expects(b = 2)
  stub_func <- stub_of_simpf$f
  
  expect_null(stub_func(1, 2, 3, c = 4))
  expect_null(stub_func(c = 'b', 2, a = 1, list(a = 1)))
  
  expect_error(stub_func(2, 3, 3, c = 4), not_expected_error)
  expect_error(stub_func(2, 3, 3), not_expected_error)
  expect_error(stub_func(c = 4, a = 3, 1, 2), not_expected_error)
  expect_error(stub_func(a = 3, 1, 2), not_expected_error)
})

test_that('expects & returns: Always checks for expected arguments and returns the specified value', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$expects(c = list(a = '1'))
  stub_of_simpf$returns('a')
  stub_func <- stub_of_simpf$f
  
  expect_equal(stub_func(1, 2, 3, c = list(a = '1')), 'a')
  expect_equal(stub_func(c = list(a = '1'), 2, a = 1, 3), 'a')
  
  expect_error(stub_func(2, 3, 3, c = 4), not_expected_error)
  expect_error(stub_func(2, 3, 3), not_expected_error)
  expect_error(stub_func(c = 4, a = 3, 1, 2), not_expected_error)
})

test_that('expects & throws: Always checks for expected arguments and throws error with specified msg', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$expects(c = 'p')
  stub_of_simpf$throws('err msg')
  stub_func <- stub_of_simpf$f
  
  expect_error(stub_func(1, 2, 3, c = 'p'), 'err msg')
  expect_error(stub_func(c = 'p', list(1, 2, 3), a = 1, 3), 'err msg')
  
  expect_error(stub_func(2, 3, 3, c = 4), not_expected_error)
  expect_error(stub_func(2, 3, 3), not_expected_error)
  expect_error(stub_func(c = 4, a = 3, 1, 2), not_expected_error)
})