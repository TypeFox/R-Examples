context("EDR")

test_that("The EDR function is calculated correctly", {
  
  #FIrst example (no window)
  q<-c(0) 
  r<-c(1, 2)
  s<-c(2, 3, 3)

  #EDR distance between x and y with no window is  (see )
  expect_equal(EDRDistance(q, r, epsilon=1), 1)
  expect_equal(EDRDistance(r, s, epsilon=1), 1)
  expect_equal(EDRDistance(q, s, epsilon=1), 3)

  
  #Second example, window size=1
  r <- c(1, 9, 1, 7, 5, 12)
  s <- c(9, 6, 3, 11, 2, 1)
  
 #With sigma=1 the cost matrix is
 # 1 2 # # # #
 # 1 2 3 # # #
 # # 3 3 4 # #
 # # # 4 4 5 # 
 # # # # 5 5 6
 # # # # # 6 6
 
 # Two minimal mappings 
 # 1-1, 2-2, 3-3, 4-4, 5-5, 6-6
 # ?-1, 1-2, 2-? 3-3, 4-4, 5-5, 6-6 
 
 expect_equal(EDRDistance(r, s, epsilon=1, sigma=1), 6)
  
 #With sigma=2 the cost matrix is
 
 # 1 2 3 # # #
 # 1 2 3 4 # #
 # 2 2 3 4 4 #
 # # 2 3 4 5 5 
 # # # 3 4 5 6
 # # # # 3 4 5
 
 #Example of minimal match
 
 #1-?, 2-1, 3-?, 4-2, 5-3, 6-4, 6-5, 6-6 

 expect_equal(EDRDistance(r, s, epsilon=1, sigma=2), 5)
})


test_that("Exceptions in EDR distance", {
  
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  epsilon<-1
  expect_equal(EDRDistance(x, y, epsilon), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(EDRDistance(x, y, epsilon), NA)
  
  x <- as.numeric(c())
  expect_equal(EDRDistance(x, y, epsilon), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(EDRDistance(x, y, epsilon), NA)
  
  x <- c(1, 2, 3, 4)
  epsilon <- "a"
  expect_equal(EDRDistance(x, y, epsilon), NA)
  
  epsilon <- -1
  expect_equal(EDRDistance(x, y, epsilon), NA)
  
  epsilon <-1 
  sigma <- -1
  expect_equal(EDRDistance(x, y, epsilon, sigma), NA)
  
  sigma <- 10
  expect_equal(EDRDistance(x, y, epsilon, sigma), NA)
  
  x <- c(1, 2, 3, 4, 5, 6)
  sigma <- 5
  expect_equal(EDRDistance(x, y, epsilon, sigma), NA)
  
  x <- c(1, 2)
  sigma <- 1
  expect_equal(EDRDistance(x, y, epsilon, sigma), NA)
  
})
