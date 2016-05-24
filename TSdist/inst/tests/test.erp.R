context("ERP")

test_that("The ERP function is calculated correctly", {
  
  #FIrst example (no window)
  q<-c(0) 
  r<-c(1, 2)
  s<-c(2, 3, 3)
  
  #erp distance between x and y with no window is  (see )
  expect_equal(ERPDistance(q, r, g=0), 3)
  expect_equal(ERPDistance(r, s, g=0), 5)
  expect_equal(ERPDistance(q, s, g=0), 8)
  
  
  #Second example
  r <- c(1, 9, 1, 7, 5, 12)
  s <- c(9, 6, 3, 11, 2, 1)
  
  #With window size 1, the cost matrix is 
  
  #  8  1  #  #  #  #
  # 14  7  6  #  #  #
  #  # 10  9 10  #  #
  #  #  # 20 15 16  # 
  #  #  #  # 17 18 26
  #  #  #  #  # 19 27
  
  #An minimal matching:
  #1-?, 2-1, 3-2, 4-3, 5-4, 6-5, ?-6
  
  expect_equal(ERPDistance(r, s, g=0, sigma=1), 27)
  

  #With window size 2, the cost matrix is
  
  #  8  1  2  #  #  #
  # 14  7  6  3  #  #
  # 17 10  9  6  5  #
  #  # 19 20 13 12  6 
  #  #  # 20 15 14  8
  #  #  #  # 16 15  9
  
  #An optimal matching:
  #1-?, 2-1, 3-?, 4-2, 5-3, 6-4, ?-5, ?-6
  
  expect_equal(ERPDistance(r, s, g=0, sigma=2), 9)
})


test_that("Exceptions in erp distance", {
  
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  g<-1
  expect_equal(ERPDistance(x, y, g), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(ERPDistance(x, y, g), NA)
  
  x <- as.numeric(c())
  expect_equal(ERPDistance(x, y, g), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(ERPDistance(x, y, g), NA)
  
  x <- c(1, 2, 3, 4)
  g <- "a"
  expect_equal(ERPDistance(x, y, g), NA)
  
  g <-1 
  sigma <- -1
  expect_equal(ERPDistance(x, y, g, sigma), NA)
  
  sigma <- 10
  expect_equal(ERPDistance(x, y, g, sigma), NA)
  
  x <- c(1, 2, 3, 4, 5, 6)
  sigma <- 5
  expect_equal(ERPDistance(x, y, g, sigma), NA)
  
  x <- c(1, 2)
  sigma <- 1
  expect_equal(ERPDistance(x, y, g, sigma), NA)
  
})