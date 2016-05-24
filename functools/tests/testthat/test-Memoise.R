library(functools)
context("Memoise()")

fib_R <- function(n) {
  if (n == 0) return(0)
  if (n == 1) return(1)
  return (fibR(n - 1) + fibR(n - 2))
}

fib_R_memoised <- Memoise(fib_R)

x <- round(runif(100, 0, 100))

runif_memoised <- Memoise(runif)

test_that("Produces the correct output.", {
  expect_equal(runif_memoised(100), runif_memoised(100))
})

test_that("Produces the correct output type.", {
  expect_is(Memoise(Identity), "function")
})

test_that("Produces the correct errors.", {
  expect_error(Memoise()())
})
