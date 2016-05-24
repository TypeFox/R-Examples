context("Testing Walkr Output")

## Tests that Walkr returns the desired 
## output

test_that("Testing Walkr Output", {
  
  set.seed(314)
  
  ## Simple 3D simplex
  
  A1 <- matrix(c(1,0,1), ncol = 3)
  b1 <- 0.5
  
  ## Random 50D A 
  
  A2 <- matrix(sample(c(1,0), 50, replace = T), ncol = 50)
  b2 <- 0.6

  ## sampling
  ## suppress warnings because walkr indicates mixing might not be enough
  
  z1_har <- suppressWarnings(walkr(A = A1, b = b1, points = 50, 
                                   method = "hit-and-run", chains = 5))
  z1_dikin <- suppressWarnings(walkr(A = A1, b = b1, points = 50, 
                                     method = "dikin", chains = 5))
  
  z2_har <- suppressWarnings(walkr(A = A2, b = b2, points = 50, 
                                   method = "hit-and-run", chains = 5))
  z2_dikin <- suppressWarnings(walkr(A = A2, b = b2, points = 50, 
                                     method = "dikin", chains = 5))
  
  ## check that we're returning a matrix
  
  expect_equal(class(z1_har), "matrix")
  expect_equal(class(z1_dikin), "matrix")
  expect_equal(class(z2_har), "matrix")
  expect_equal(class(z2_dikin), "matrix")
  
  ## check that we the right dimensions
  
  expect_equal(dim(z1_har), c(3, 50))
  expect_equal(dim(z1_dikin), c(3, 50))
  expect_equal(dim(z2_har), c(50, 50))
  expect_equal(dim(z2_dikin), c(50, 50))
 
  ## check that we have no NA's in the output
  
  expect_true(!any(is.na(z1_har)))
  expect_true(!any(is.na(z1_dikin)))
  expect_true(!any(is.na(z2_har)))
  expect_true(!any(is.na(z2_dikin)))
  
  ## higher dimensions should spit out warning that mixing 
  ## is not good
  
  expect_warning(walkr(A = A2, b = b2, points = 50, 
                                   method = "hit-and-run", chains = 5), 
                "there are parameters with rhat > 1.1, you may want to run your chains for longer")
  
  
  
  
})
