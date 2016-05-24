# testthat.R -- automated test cases 
#
# MS 
# Created: Feb 8, 2015

require("testthat")
require("Ckmeans.1d.dp")

test_that("Ckmeans.1d.dp() finding the number of clusters", {

  x <- c(.9, 1, 1.1, 1.9, 2, 2.1)
  result <- Ckmeans.1d.dp(x)
  expect_equal(result$size, c(3,3))

  x <- rev(x)
  result <- Ckmeans.1d.dp(x)
  expect_equal(result$size, c(3,3))
  
  x <- rep(1, 100)
  result <- Ckmeans.1d.dp(x)
  expect_equal(result$size, 100)
  
  x <- 1:10
  result <- Ckmeans.1d.dp(x, k=c(1,10))
  expect_equal(result$size, 10)
    
})

test_that("Ckmeans.1d.dp() given the number of clusters", {
  
  x <- c(-1, 2, -1, 2, 4, 5, 6, -1, 2, -1)
  result <- Ckmeans.1d.dp(x, 3)
  expect_equal(result$size, c(4,3,3))
  expect_equal(result$cluster, c(1,2,1,2,3,3,3,1,2,1))
  expect_equal(result$centers, c(-1, 2, 5))
  expect_equal(result$withinss, c(0,0,2))

})