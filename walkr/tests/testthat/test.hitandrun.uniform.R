context("Testing hit-and-run uniformity")

test_that("Testing hit-and-run uniformity", {
  
  set.seed(314)
  ## Simplest possible case 
  
  A <- matrix(1, ncol = 3)
  b <- 1
  
  ## In this case, we are sampling from the 3D simplex
  
  ## Draw a sample of 1,000. If the sample is truly random,
  ## then each variable should have a 50/50 chance of being above 0.5. It 
  ## follows that the number of times, out of a thousand, that it is above 0.5 
  ## is provided by a binomial distrbution. So, we check that number of times 
  ## above 0.5, for both variables, is less than the 99th percentile of the
  ## binomial distribution.
  
  z <- walkr(A = A, b = b, points = 1000, method = "hit-and-run")
  z1 <- z[1,]
  z2 <- z[2,]
  standard.dev <- sd(z1-z2)
  
  conf <- qnorm(p = c(0.01, 0.99), mean = 0, sd = standard.dev / sqrt(1000))
  
  actual <- mean(z1-z2)
  
  expect_true(actual > conf[1])
  expect_true(actual < conf[2])
 
  ### A 5D Simplex
  
  A <- matrix(1, ncol = 5)
  b <- 1 
  
  ## the thinning is needed
  
  z <- walkr(A = A, b = b, points = 10000, thin = 100, method = "hit-and-run")
  
  ## should expect that the mean of x_1, x_2 be roughly the same as the mean of x_4, x_5
  
  s1 <- z[1,] + z[2,]
  s2 <- z[4,] + z[5,]
  standard.dev <- sd(s1-s2)
  
  ## expected value of their difference is zero
  ## normal approximation applicable here
  ## 99% confidence interval
  
  conf_interval.99 <- qnorm(p = c(0.01, 0.99), mean = 0, sd = standard.dev / sqrt(10000))
  
  expect_true(mean(s1-s2) >= conf_interval.99[1])
  expect_true(mean(s1-s2) <= conf_interval.99[2])
   
})