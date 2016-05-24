context("Cross Correlation")

test_that("Cross correlation function is calculated correctly", {

  x <- c(1:5)
  y <- c(2,3,1,4,5)
  
  #Mean and standard deviation of x and y
  meanx <- mean(x)
  meany <- mean(y)
  sdx <- sqrt(sum((x - meanx)^2))
  sdy <- sqrt(sum((y - meany)^2))
  
  
  #We manually calculate the cross correlation of these two
  #series until lag 2
  cc0 <- sum((x - meanx) * (y - meany)) / (sdx * sdy)
  cc1 <- sum((x[1:4] - meanx) * (y[2:5] - meany)) / (sdx * sdy)
  cc2 <- sum((x[1:3] - meanx) * (y[3:5] - meany)) / (sdx * sdy)
  
  #We calculate the distance manually
  d <- sqrt((1 - cc0 ^ 2) / (cc1 ^ 2 + cc2 ^ 2))
  
  expect_equal(CCorDistance(x, y, lag.max=2), d)
})

test_that("Exceptions in cross correlation distance", {
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  expect_equal(CCorDistance(x, y), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(CCorDistance(x, y), NA)
  
  x <- as.numeric(c())
  expect_equal(CCorDistance(x, y), NA)
  
  
  x <- c(1, 2, NA, 3)
  expect_equal(CCorDistance(x, y), NA)
  
  x <- c(1, 2, 3, 4)
  lag.max <- -1
  expect_equal(CCorDistance(x, y, lag.max=lag.max), NA)
  
  lag.max <- 10
  expect_equal(CCorDistance(x, y, lag.max=lag.max), NA)
  
  x <- c(1, 2, 3, 4, 5, 6)
  lag.max <- 5
  expect_equal(CCorDistance(x, y, lag.max=lag.max), NA)
  
})
