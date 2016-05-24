context("Dissim")

test_that("The Dissim function is calculated correctly", {
  
  #FIrst example
  x <- c(5, 5, 5, 5, 5, 5, 5)
  y <- c(2, 2, 2, 2, 2, 2, 2)
  tx <- c(1:7)
  ty <- c(1:7)
  
  # Euclidean distance between x and y is always 3 and 
  # the integral of 3 between 1 and 7 is 3*6=18
  expect_equal(DissimDistance(x, y, tx, ty), 18)
  
  
  x <- c(0, 0, 0, 0, 0)
  y <- seq(0, 1, 1/4)
  tx <- seq(0, 1, length.out=length(x))
  ty <- seq(0, 1, length.out=length(y))
   
  # Euclidean distance between x and y in time Euc(t)=sqrt(t^2)
  # and if t is always positive Euc(t)=t.
  
  # The integral of this is t^2/2, which in [0,1] takes a value of
  # 1^2/2-0^2/2=0.5
  expect_equal(DissimDistance(x, y, tx, ty), 0.5)
  
  # For two series of different length (value of distance is the same as in the 
  # previous case)
  y <- y[-3]
  ty <- ty[-3]
  expect_equal(DissimDistance(x, y, tx, ty), 0.5)
  
  # For two series of same length but different time indices (also based on 
  # previous example)
  x <- x[-2]
  tx <- tx[-2]
  expect_equal(DissimDistance(x, y, tx, ty), 0.5)
  
  
  # If we don't define tx and ty, they are defined by dividing [0,1] into an
  # equally spaced grid
  expect_equal(DissimDistance(x, y), 0.5)
  expect_equal(DissimDistance(x, y, tx=tx), 0.5)
  expect_equal(DissimDistance(x, y, ty=ty), 0.5)
})


test_that("Exceptions in Dissim distance", {
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  n <- 2
  expect_equal(DissimDistance(x, y), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(DissimDistance(x, y), NA)
  
  x <- as.numeric(c())
  expect_equal(DissimDistance(x, y), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(DissimDistance(x, y), NA)
  
  x <- c(1, 2, 3, 4)
  tx <- c(1, 2, 3, 4)
  ty <- c(2, 3, 3.5, 4)
  expect_equal(DissimDistance(x, y, tx, ty), NA)
  
  tx <- c(1, -2, 3, 4)
  ty <- c(1, 2, 3, 4)
  expect_equal(DissimDistance(x, y, tx, ty), NA)
  
  tx <- c(1, 3, 2, 4)
  expect_equal(DissimDistance(x, y, tx, ty), NA)
  
  tx <- c(1, 2, 4)
  expect_equal(DissimDistance(x, y, tx, ty), NA)
})
