context("TSDistances")

test_that("The TSDistances function is checked", {
  
  x <- cumsum(rnorm(100))
  y <- x
  
  # If x and y are given as ts objects
  
  x1 <- as.ts(x)
  x2 <- as.ts(y)
  expect_equal(TSDistances(x1, x2, "euclidean"), 0)


  # If x and y are given as zoo objects
  
  x1 <- as.zoo(x)
  x2 <- as.zoo(y)
  expect_equal(TSDistances(x1, x2, "euclidean"), 0)
  
  # If x and y are given as xts objects
  data(zoo.series1)
  x1 <- as.xts(zoo.series1) 
  x2 <- as.xts(zoo.series1) 
  expect_equal(TSDistances(x1, x2, "euclidean"), 0)
  
  
})
