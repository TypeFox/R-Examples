library(testthat)

simpf <- function(a = 1, b, d = paste, ...) return(5)
not_expected_error <- 'Function is called with arguments different from expected!'

test_that('Returns the specified value if expected arguments are part of the function call', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withArgs(b = 2)$returns(10)
  stub_func <- stub_of_simpf$f

  expect_equal(stub_func(1, 2, 3, f = 4), 10)
  expect_equal(stub_func(3, 2, 5), 10)
  expect_null(stub_func(3, 4, 5))
})

test_that('Returns the specified value if expected arguments are part of the function call', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withArgs(f = 2)$returns(10)
  stub_func <- stub_of_simpf$f
  
  expect_equal(stub_func(1, 2, 3, f = 2), 10)
  expect_null(stub_func(1, 2, 3, g = 2))
})


test_that('Throws error with specified message if called with the exact arguments specified', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withArgs(b = 2)$throws('pqrs')
  stub_func <- stub_of_simpf$f

  expect_error(stub_func(1, 2, 3, f = 4), 'pqrs')
  expect_error(stub_func(3, 2, 5), 'pqrs')
  expect_null(stub_func(3, 4, 5))
})

test_that('Throws error with specified message if called with the exact arguments specified', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withArgs(f = 2)$throws('pqrs')
  stub_func <- stub_of_simpf$f
  
  expect_error(stub_func(1, 2, 3, f = 2), 'pqrs')
  expect_null(stub_func(3, 4, 5))
  expect_null(stub_func(1, 2, 3, g = 2))
})

test_that('It does the right thing even when there are multiple expectations - withArgs.return/throw and default return', {
  stub_of_simpf <- stub(simpf)
  stub_of_simpf$withArgs(b = 2)$returns(10)
  stub_of_simpf$withArgs(c = 5, d = 8)$throws('xyz')
  stub_of_simpf$withArgs(b = 3)$returns(20)
  stub_of_simpf$throws('pwrs')
  stub_func <- stub_of_simpf$f

  expect_equal(stub_func(1, 2, 3, c = 4), 10)
  expect_error(stub_func(1, 1, 3, c = 5), 'pwrs')
  expect_error(stub_func(1, 1, 3, c = 5, d = 8), 'xyz')
  expect_equal(stub_func(1, 3, 4, c = 4), 20)
  expect_error(stub_func(9, 9, 9, c = 9), 'pwrs')
})

