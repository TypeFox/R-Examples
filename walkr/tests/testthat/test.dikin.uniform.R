context("Testing dikin uniformity")

test_that("Testing dikin uniformity", {
  
  ## The Simple 3D simplex
  set.seed(314)
  A <- matrix(1, ncol = 3)
  b <- 1
  
  ## We should expect to see that all coordinates have
  ## the same value. Note that the expected value of each of the coordinates 
  ## is 1/3 
  
  ## construct confidence interval
  ## standard error divided by sample size
  
  z <- walkr(A = A, b = b, points = 1000, method = "dikin")
  
  conf1 <- qnorm(p = c(0.01, 0.99), mean = 1/3, sd = sqrt(sd(z[1,])/1000))
  conf2 <- qnorm(p = c(0.01, 0.99), mean = 1/3, sd = sqrt(sd(z[2,])/1000))
  conf3 <- qnorm(p = c(0.01, 0.99), mean = 1/3, sd = sqrt(sd(z[3,])/1000))
  
  
  ## all should fall within confidence interval
  
  expect_true(mean(z[1,]) <= conf1[2])
  expect_true(mean(z[1,]) >= conf1[1])
  expect_true(mean(z[2,]) <= conf2[2])
  expect_true(mean(z[2,]) >= conf2[1])
  expect_true(mean(z[3,]) <= conf3[2])
  expect_true(mean(z[3,]) >= conf3[1])
  
  
  
})