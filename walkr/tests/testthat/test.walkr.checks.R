context("Testing Walkr")

## Tests that Walkr throws the correct error messages when 
## input is incorrect

test_that("Testing Walkr Checks and Input", {
  
  set.seed(314)
  
  ## check that A is a matrix
  
  A <- data.frame(1, ncol = 3)
  b <- 1
  #method doesn't matter,
  expect_error( walkr(A = A, b = b, points = 10), "A needs to be a matrix")
  
  ## check that A is underdetermined
  
  A <- matrix(stats::runif(9,0,1), ncol = 3)
  b <- 1
  expect_error(walkr(A = A, b = b, points = 10), "A must be underdetermined")
  
  ## check that A and b matches
  
  A <- matrix(1, ncol = 3)
  b <- c(1, 0.1)
  expect_error(walkr(A = A, b = b, points = 10), "Dimensions of A and b don't match")
  
  ## A and b need to be numeric only
  
  A <- matrix(1, ncol = 3)
  b <- "hello world"
  expect_error(walkr(A = A, b = b, points = 10), "b needs to be a numeric vector")
  
  A <- matrix("hello", ncol = 3)
  b <- 1
  expect_error(walkr(A = A, b = b, points = 10), "A needs to contain numbers only")
  
  ## Ill-defined method
  
  A <- matrix(1, ncol = 3)
  b <- 1
  expect_error(walkr(A = A, b = b, points = 10, method = "ball-search"), 
               "Method must be hit-and-run or dikin")
  
  ## Correct points, thin, chains parameter
  A <- matrix(1, ncol = 3)
  b <- 1
  expect_error(walkr(A = A, b = b, points = 0), 
               "points, thin, chains must be geq 1")
  expect_error(walkr(A = A, b = b, points = 1, thin = 0), 
               "points, thin, chains must be geq 1")
  expect_error(walkr(A = A, b = b, points = 1, thin = 1, chains = 0), 
               "points, thin, chains must be geq 1")
  
  ## burn non-negative
  expect_error(walkr(A = A, b = b, points = 1, burn = -1), 
               "burn must be non-negative")
 
  
  
})
