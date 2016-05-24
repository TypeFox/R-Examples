context("STS")

test_that("The STS function is calculated correctly", {
  
  x <- c(1, 1, 1, 1, 1)
  y <- c(2, 2, 2, 2, 2)
  
  # The function calculates the sum of the squared differences 
  # between the derivatives. In this case, the derivatives are 
  # the same throughout the series so the sts distance should be 0
  
  expect_equal(STSDistance(x, y), 0)
  
  x <- c(0, 1, 2, 3, 4, 5)
  y <- c(5, 4, 3, 2, 1, 0)
  
  # The derivative of the first series is always 1.
  # The derivative of the first series is always -1.
  # The distance is thus 
  #   sqrt((1-(-1))^2*(number of points-1)=sqrt(4*5)=sqrt(20)

  expect_equal(STSDistance(x, y), sqrt(20))
  
  # If we define one of the indices but not the other, they will be set to 
  # equal.
  expect_equal(STSDistance(x, y, tx=c(1:length(x))), sqrt(20))
  expect_equal(STSDistance(x, y, ty=c(1:length(y))), sqrt(20))
  
})


test_that("Exceptions in sts distance", {
  
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  expect_equal(STSDistance(x, y), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(STSDistance(x, y), NA)
  
  x <- c(1, 2)
  expect_equal(STSDistance(x, y), NA)
  
  x <- as.numeric(c())
  expect_equal(STSDistance(x, y), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(STSDistance(x, y), NA)
  
  x <- c(1, 2, 3, 4)
  tx <- c(1, -1, 2, 3)
  ty <- c(1, 2, 3, 4)
  expect_equal(STSDistance(x, y, tx, ty), NA)
  
  tx <- c(1, 3, 4, 5)
  expect_equal(STSDistance(x, y, tx, ty), NA)
  
  tx <- c(1, 3, 2, 4)
  ty <- c(1, 3, 2, 4)
  expect_equal(STSDistance(x, y, tx, ty), NA)
  
  tx <- c(1, 2, 3)
  ty <- c(1, 2, 3, 4)
  expect_equal(STSDistance(x, y, tx, ty), NA)
  
})