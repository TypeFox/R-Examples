# 
# 


context('Testing inahull_cpp')

test_that("inahull_cpp works", { 
  library(rollply)
  library(alphahull)
  
  f = 1
  xy <- as.data.frame(matrix(rnorm(1000)*f, ncol = 2))
  hull <- alphahull::ahull(xy, alpha = .4 * f)

  xy.test <- expand.grid(x = seq(-1, 1, length.out = 10)*f*2, 
                         y = seq(-1, 1, length.out = 10)*f*2)

  test <- apply(xy.test, 1, function(X) alphahull::inahull(hull, X))
  test2 <- inahull_cpp_multiple(hull, xy.test[[1]], xy.test[[2]])
  expect_true( all.equal(test, test2) )
  
})
