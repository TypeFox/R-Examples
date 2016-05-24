#require(testthat); library("gapfill", lib.loc = "../../../lib/")
context("test-Score")

test_that("Score_R",{
  Score_R <- function(mat){
      out <- rep(NA, dim(mat)[2])
      for(i in 1:dim(mat)[2]){
          others <- mat[,-i,drop=FALSE] 
          ref <- array(mat[,i], dim(others)) 
          diff <- ref - others
          out[i] <- mean(apply(diff > 0, 2, mean,
                               na.rm=TRUE),
                         na.rm=TRUE)
      }
      out
  }

  x <- array(c(
1.32, 0.25, -0.76, -1.07, -0.04, 0.96, -0.79, 1.26, 0.2, -0.58, 
-1.16, -0.6, 0.47, 1.1, -1.95, 0.33, -1.01, 1.27, 0.43, -1.92, 
0.48, 0.26, 0.22, 0.06, 0.73, -0.18, 0.78, 0.38, 0.41, 0.12, 
0.4, 0.22, 0.45, 0.18, -0.2, 1.67, -0.43, -0.5, -0.83, 0.83, 
-1.04, 0.82, 0.8, -1.02, -1.33, -0.11, -0.81, -0.23, -0.92, -0.66
), c(10, 5))
  expect_equal(Score(x), Score_R(x))

  x[c(3L, 4L, 9L, 10L, 15L, 17L, 18L, 23L, 24L, 31L, 33L, 35L, 39L, 
      43L, 49L)] <- NA
  expect_equal(Score(x), Score_R(x))

  x[1,] <- NA
  expect_equal(Score(x), Score_R(x))

  x[,1] <- NA
  expect_equal(Score(x), Score_R(x))
 
})

test_that("Score_errors",{
  expect_error(Score(1), "not a matrix")
  expect_equal(Score(matrix(1)), NaN)
  expect_equal(Score(matrix(1L)), NaN)
  expect_equal(Score(matrix(c(TRUE,TRUE,FALSE,TRUE), 2)), c(.5, 0))
  expect_equal(Score(matrix(c(1,1,1,1,0,1), 2)), c(.25,.25,0))
})

