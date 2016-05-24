context("Test get_rates")

test_that("test for get_rates",{
  ## US treasury yield curve from April 15, 2014
   truth.1 <- c(0.001517, 0.001923, 0.002287, 0.003227, 0.005465, 0.005105, 0.009265,
                0.013470, 0.017150, 0.020160, 0.022630, 0.024580, 0.026265, 0.027590,
                0.029715, 0.031820, 0.033635, 0.034420, 0.034780)
  
  #rates column from the package's getRates function
  
  rates <- get_rates(date = as.Date("2014-04-15"), currency = "USD")$rate
 
  stopifnot(all.equal(truth.1, rates))
})
