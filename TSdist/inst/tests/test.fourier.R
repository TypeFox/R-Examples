context("Fourier")

test_that("The Fourier function is calculated correctly", {

  #We use the linearity condition of fourier transforms
  x <- rnorm(10)
  y <- x * 5
    
  #We calculate the Fourier Transforms
  fx <- fft(x)
  fy <- fft(y)
  
  #We know that fx*5=fy, we can calculate the Euclidean distance
  #betwee fx and fy in terms of fx:
  real <- Re(fx)
  imaginary <- Im(fx)
  
  #Only considering the first 2 coefficients
  d <- 4 * sqrt(sum(real[1:2] ^ 2 + imaginary[1:2] ^ 2))
  expect_equal(FourierDistance(x, y, n=2), d)
  
  #Only considering the first 4 coefficients
  d <- 4 * sqrt(sum(real[1:4] ^ 2 + imaginary[1:4] ^ 2))
  expect_equal(FourierDistance(x, y, n=4), d)

  
  })


test_that("Exceptions in Fourier distance", {
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  n <- 2
  expect_equal(FourierDistance(x, y, n), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(FourierDistance(x, y, n), NA)
  
  x <- c(1)
  expect_equal(FourierDistance(x, y, n), NA)
  
  x <- c(1, 2)
  expect_equal(FourierDistance(x, y, n), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(FourierDistance(x, y, n), NA)
  
  x <- c(1, 2, 3, 4)
  n <- 5
  expect_equal(FourierDistance(x, y, n), NA)
  
  x <- c(1, 2, 3, 4)
  n <- 4
  expect_warning(FourierDistance(x, y, n))
    
})
