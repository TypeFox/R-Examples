context("test-Array2Matrix.R")
test_that("IndexTwoOne",{
  IndexOneTwo_alt <- function(index, dim){
      tmp <- array(as.numeric(1:prod(dim)), dim)
      c(which(index == tmp, arr.ind = TRUE))
  }
 
  for(i in 1:2){
      for(j in 1:3){
         expect_equal(IndexOneTwo_alt(IndexTwoOne(c(i,j), c(2,3)),
                                       c(2,3)), c(i,j))
         expect_equal(IndexOneTwo_alt(IndexTwoOne(c(j,i), c(3,2)),
                                       c(3,2)), c(j,i))
     }
  }
})

test_that("IndexOneFour",{
              IndexOneFour_alt <- function(index, dim){
                  tmp <- array(1:prod(dim), dim)
                  stopifnot(length(index) == 1)
                  c(which(index == tmp, arr.ind = TRUE))
              }

              dim <- rep(1,4)
              for(i in 1:prod(dim)){
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
              }
              dim <- c(3,3,3,3)
              for(i in 1:prod(dim)){
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
              }
              dim <- 2:5
              for(i in 1:prod(dim)){
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
              }
              dim <- 5:2
              for(i in 1:prod(dim)){
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
                  expect_equal(IndexOneFour_alt(i, dim),
                               IndexOneFour(i, dim))
              }
})


test_that("Array2Matrix",{
   expect_error(Array2Matrix(array(1:12, c(1,3,4))))
   
   x <- array(1:12, c(1,3,4,1))
   expect_equal(Array2Matrix(x), array(1:12, c(3,4)))

   attr(x, "mp") <- c(1,2,2,1)
   tmp <- Array2Matrix(x)
   attr(tmp, "mp") <- NULL
   expect_equal(tmp, array(1:12, c(3,4)))

   y <- array(1:12, c(3,4))
   attr(y, "mp") <- c(2, 2)
   expect_equal(y, Array2Matrix(x))
})
