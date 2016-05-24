context("Keogh_LB")

test_that("The Keogh_LB function is calculated correctly", {
  
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 3, 5, 2, 4)
  
  #With a window of size 1, the upper and lower envelope is
  #the series itself
  u <- x
  l <- x
  
  #Now, y is stricly larger than u in points 2 (3>2) and 3 (5>3)
  case1 <- which(y>u)
  
  #For these points, we calculate the contribution to the 
  #Keogh_LB distance as: (y-u)^2
  
  d1 <- (y[case1] - u[case1]) ^ 2
  
  #y is strictly lower than l in points 4 (2<4)  and 5 (4<5).
  case2 <- which(y < l)
  
  #For these points, we calculate the contribution to the 
  #Keogh_LB distance as: (y-l)^2
  d2 <- (y[case1] - l[case1]) ^ 2
  
  #For the points that lie in between u and l, the contribution
  #to the lower bound is 0. 
  
  #We sum all contributions and take its root to obtain the lower bound:
  d<-sqrt(sum(c(d1, d2)))
  
  #lb.keogh distance between x and y with no window is  (see )
  expect_equal(LBKeoghDistance(x, y, window.size=1), d)

  })


test_that("Exceptions in lb.keogh distance", {
  
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  expect_equal(LBKeoghDistance(x, y, window.size=1), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(LBKeoghDistance(x, y, window.size=1), NA)
  
  x <- c(1, 2)
  expect_equal(LBKeoghDistance(x, y, window.size=1), NA)
  
  x <- as.numeric(c())
  expect_equal(LBKeoghDistance(x, y, window.size=1), NA)
  
  x <- c(1, 2, NA, 3)
  expect_equal(LBKeoghDistance(x, y, window.size=1), NA)
  
  x <- c(1, 2, 3, 4)
  expect_equal(LBKeoghDistance(x, y, window.size=2), NA)
  
  expect_equal(LBKeoghDistance(x, y, window.size=11), NA)
  
})