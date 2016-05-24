context("Testing Affine Transformation")

test_that("Affine Transformation", {
  
  set.seed(314)
  
  ## a random 20-D polytope
  
  n <- 20
  A <- matrix(rep(1,n), ncol = n, nrow = 1)
  A <- rbind(A, sample(c(1,0), n, replace = T))
  A <- rbind(A, sample(c(1,0), n, replace = T))
  A <- rbind(A, sample(c(1,0), n, replace = T))
  b <- c(1, 0.7, 0.2, 0.05)
  
  ## we want to check that our solution basis is indeed correct
  
  z <- walkr:::complete_solution(A, b)
  
  ## since any linear combination of the null space basis span the null space
  ## thus, particular + (homogeneous %*% any linear combination) = b 
  
  expect_true(all((A %*% (z[[1]] + 575 * z[[2]][ ,1] + 300 * z[[2]][ ,2])) - b <= 1e-5))
  
  
})