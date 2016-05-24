context('pp')

test_that('it can interpolate a simple example', {
  x <- 1
  expect_equal(pp("1 + #{x}"), "1 + 1")
})

test_that('it can interpolate a non-name expression', {
  x <- 1
  expect_equal(pp("1 + #{x + 1}"), "1 + 2")
})

