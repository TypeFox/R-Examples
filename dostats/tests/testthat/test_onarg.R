test_that("testing onarg", {
  expect_is(onarg(runif,'max'), 'function')
  expect_true(all(onarg(runif, 'max')(1:10,10) < 1:10))
  expect_equal(length(letters %contains% 'a'), 1)
  expect_is(letters %contains% 'a', 'logical')
})
