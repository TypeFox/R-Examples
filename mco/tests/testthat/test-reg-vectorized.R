##
## Make sure vectorized return values are handled correctly
##

context("reg-vectorized")

f <- function(x)
  c(1, 2, 3)

vf <- function(X) {
  ## Number of different parameter combinations passed in:
  n <- nrow(X)
  matrix(1:3, nrow=3, ncol=n)
}

test_that("Auto-vectorization works", {
  r <- nsga2(f, idim=2, odim=3,
             lower.bounds=c(-1, -1), upper.bounds=c(1, 1),
             popsize=12, generations=1)
  objectives <- r$value
  ## Note this is transposed compared to the expected return value of f!
  expected_value <- matrix(as.numeric(1:3), nrow=12, ncol=3, byrow=TRUE)
  expect_equal(objectives, expected_value)
})

test_that("Vectorization works", {
  r <- nsga2(vf, idim=2, odim=3, vectorized=TRUE,
             lower.bounds=c(-1, -1), upper.bounds=c(1, 1),
             popsize=12, generations=1)
  objectives <- r$value
  ## Note this is transposed compared to the expected return value of f!
  expected_value <- matrix(as.numeric(1:3), nrow=12, ncol=3, byrow=TRUE)
  expect_equal(objectives, expected_value)
})
