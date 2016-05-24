context("Lp")

test_that("Euclidean function is calculated correctly", {
  x <- c(1, 2, 3, 4)
  y <- c(3, 4, 1, 2)
  #The euclidean distance of these two functions is sqrt(2^2+2^2+2^2+2^2)=sqrt(16)=4
  expect_equal(EuclideanDistance(x, y), 4)
})

test_that("Minkowski function is calculated correctly", {
  x <- c(1, 2, 3, 4)
  y <- c(3, 4, 1, 2)
  #The Minkowski distance of these two functions is (2^3+2^3+2^3+2^3)^1/3=(32)^1/3
  expect_equal(MinkowskiDistance(x, y, 3), 32^(1/3))
})

test_that("Infinite norm function is calculated correctly", {
  x <- c(1, 2, 3, 4)
  y <- c(3, 4, 1, 2)
  #The euclidean distance of these two functions is max(2, 2, 2, 2)=2
  expect_equal(InfNormDistance(x, y), 2)
})

test_that("Manhattan function is calculated correctly", {
  x <- c(1, 2, 3, 4)
  y <- c(3, 4, 1, 2)
  #The euclidean distance of these two functions is 2+2+2+2=8
  expect_equal(ManhattanDistance(x, y), 8)
})
  
test_that("Euclidean. Manhattan, Minkowski and Infinite Norm distances are 
          calculated correctly using LPDistance function", {
  x <- c(1, 2, 3, 4)
  y <- c(3, 4, 1, 2)
  expect_equal(LPDistance(x, y, "euclidean"), 4)
  expect_equal(LPDistance(x, y, "minkowski", p=3), 32^(1/3))
  expect_equal(LPDistance(x, y, "manhattan"), 8)
  expect_equal(LPDistance(x, y, "infnorm"), 2)
})

test_that("Exceptions in Lp distances", {
  x <- c("a", "b", "c", "d")
  y <- c(3, 4, 1, 2)
  expect_equal(EuclideanDistance(x, y), NA)
  expect_equal(MinkowskiDistance(x, y, p=3), NA)
  expect_equal(ManhattanDistance(x, y), NA)
  expect_equal(InfNormDistance(x, y), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(EuclideanDistance(x, y), NA)
  expect_equal(MinkowskiDistance(x, y, p=3), NA)
  expect_equal(ManhattanDistance(x, y), NA)
  expect_equal(InfNormDistance(x, y), NA)
  
  x <- as.numeric(c())
  expect_equal(EuclideanDistance(x, y), NA)
  expect_equal(MinkowskiDistance(x, y, p=3), NA)
  expect_equal(ManhattanDistance(x, y), NA)
  expect_equal(InfNormDistance(x, y), NA)
  
  x <- c(1, 2)
  expect_equal(EuclideanDistance(x, y), NA)
  expect_equal(MinkowskiDistance(x, y, p=3), NA)
  expect_equal(ManhattanDistance(x, y), NA)
  expect_equal(InfNormDistance(x, y), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(EuclideanDistance(x, y), NA)
  expect_equal(MinkowskiDistance(x, y, p=3), NA)
  expect_equal(ManhattanDistance(x, y), NA)
  expect_equal(InfNormDistance(x, y), NA)
  
  x <- c(1, 2, 2, 3)
  p <- 1.2
  expect_equal(MinkowskiDistance(x, y, p), NA)
  
  p <- -1
  expect_equal(MinkowskiDistance(x, y, p), NA)
  
})
