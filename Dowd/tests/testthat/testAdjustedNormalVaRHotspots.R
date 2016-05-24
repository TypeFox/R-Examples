test_that("Adjusted Normal VaR Hotspots.",{

  # Success - 1
  vc.matrix <- matrix(c(2.5, 3.4, -1.9, 4.3, 2.3, -3.1, 4.3, -1.2, 1), 3, 3)
  mu <- c(.4, -.3, .1)
  skew <- .5
  kurtosis <- 1.2
  positions <- c(5,2,6)
  cl <- .95
  hp <- 280
  val <- AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp)
  expect_equal(c(4870.2, 1056.3, 1280.1), val, tolerance=0.001)
  
  # Success - 2
  vc.matrix <- matrix(c(2.6, 3.4, -1.9, 4.3, 2.3, -3.1, 4.3, -1.2, 
                       1, 3.4, -2, -5, -1.2, 3.2, 0, -1.2), 4, 4)
  mu <- c(-.5, -.3, -1.2, 0)
  skew <- -.4
  kurtosis <- 2.2
  positions <- c(4,1,10,3)
  cl <- .99
  hp <- 50
  val <- AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp)
  expect_equal(c(1200.9, 1738.4, -4585.4, -1973.5), val, tolerance=1)
  
  # Error - 1
  vc.matrix <- matrix(c(2.5, 3.4, -1.9, 4.3, 2.3, -3.1, 4.3, -1.2, 1), 3, 3)
  mu <- c(.4, -.3, .1)
  skew <- .5
  kurtosis <- 1.2
  positions <- c(5,2,6)
  cl <- .95
  hp <- -10
  expect_error(val <- AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp))
  
  # Error - 2
  vc.matrix <- matrix(c(2.5, 3.4, -1.9, 4.3, 2.3, -3.1, 4.3, -1.2, 1), 3, 3)
  mu <- c(.4, -.3, .1)
  skew <- .5
  kurtosis <- 1.2
  positions <- c(5,2,6)
  cl <- 1.2
  hp <- 280
  expect_error(val <- AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp))
  
  # Error - 3
  vc.matrix <- matrix(c(2.5, 3.4, -1.9, 4.3, 2.3, -3.1, 4.3, -1.2, 1), 3, 3)
  mu <- c(.4, -.3, .1)
  skew <- .5
  kurtosis <- 1.2
  positions <- c(5,2,6)
  cl <- -.95
  hp <- 280
  expect_error(val <- AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp))
  
  # Error - 4
  vc.matrix <- matrix(c(2.5, 3.4, -1.9, 4.3, 2.3, -3.1, 4.3, -1.2, 1), 3, 3)
  mu <- c(.4, -.3, .1, 1.2)
  skew <- .5
  kurtosis <- 1.2
  positions <- c(5, 2, 6)
  cl <- -.95
  hp <- 280
  expect_error(val <- AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp))
  
  # Error - 5
  vc.matrix <- matrix(c(2.5, 3.4, -1.9, 4.3, 2.3, -3.1, 4.3, -1.2, 1), 3, 3)
  mu <- c(.4, -.3, .1)
  skew <- .5
  kurtosis <- 1.2
  positions <- c(5,2,6,3)
  cl <- -.95
  hp <- 280
  expect_error(val <- AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp))
  
})