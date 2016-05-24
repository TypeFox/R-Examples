context("LCSS")

test_that("The LCSS function is calculated correctly", {
  
  #FIrst example (no window)
  q <- c(0) 
  r <- c(1, 2)
  s <- c(2, 3, 3)
  
  #lcss distance between x and y with no window is  (see )
  expect_equal(LCSSDistance(q, r, epsilon=1), 1)
  expect_equal(LCSSDistance(r, s, epsilon=1), 2)
  expect_equal(LCSSDistance(q, s, epsilon=1), 0)
  
  
  #Second example
  r <- c(1, 9, 1, 7, 5, 12)
  s <- c(9, 6, 3, 11, 2, 1)
  
  #With window size 1, only the elements 2-1 can be matched
  expect_equal(LCSSDistance(r, s, epsilon=1, sigma=1), 1)
  

  #With window size 2, elements 2-1, 4-2 and 6-4 can be matched
  
  expect_equal(LCSSDistance(r, s, epsilon=1, sigma=2), 3)
})


test_that("Exceptions in lcss distance", {
  
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  epsilon<-1
  expect_equal(LCSSDistance(x, y, epsilon), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(LCSSDistance(x, y, epsilon), NA)
  
  x <- as.numeric(c())
  expect_equal(LCSSDistance(x, y, epsilon), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(LCSSDistance(x, y, epsilon), NA)
  
  x <- c(1, 2, 3, 4)
  epsilon <- "a"
  expect_equal(LCSSDistance(x, y, epsilon), NA)
  
  epsilon <- -1
  expect_equal(LCSSDistance(x, y, epsilon), NA)
  
  epsilon <-1 
  sigma <- -1
  expect_equal(LCSSDistance(x, y, epsilon, sigma), NA)
  
  sigma <- 10
  expect_equal(LCSSDistance(x, y, epsilon, sigma), NA)
  
  x <- c(1, 2, 3, 4, 5, 6)
  sigma <- 5
  expect_equal(LCSSDistance(x, y, epsilon, sigma), NA)
  
  x <- c(1, 2)
  sigma <- 1
  expect_equal(LCSSDistance(x, y, epsilon, sigma), NA)
  
})