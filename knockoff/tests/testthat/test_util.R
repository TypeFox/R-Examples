test_that('restoring the random seed works', {
  if (!exists('.Random.seed')) RNGkind()
  seed <- .Random.seed
  first <- with_seed(0, rnorm(10))
  second <- with_seed(0, rnorm(10))
  expect_equal(.Random.seed, seed)
  expect_equal(first, second)
})