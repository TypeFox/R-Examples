library(DNAprofiles)
context("Dist checks")

test_that(desc = "Dist checks",{  
  ## die example
  dist = list(x=1:6,fx=rep(1/6,6))
  expect_identical(check.dist(dist),dist)
  dist$x[1]=0
  expect_identical(check.dist(dist),dist)
  dist$x[length(dist$x)]=Inf
  expect_identical(check.dist(dist),dist)
  dist$x[1]=-1
  expect_error(check.dist(dist))
  
  # some more fail checks
  dist1 <- list(x=1:10,fx=rep(1/10,10))
  dist1$x[5]=NA
  expect_error(check.dist(dist1))
  dist1$x[5]=NaN
  expect_error(check.dist(dist1))
  dist1$x[5]=4
  expect_error(check.dist(dist1))
    
  # unequal lengths of x and fx
  dist2 <- list(x=1:2,fx=c(1/2,1/2,0))
  expect_error(check.dist(dist2))
  
  # single event must be positive and finite  
  expect_error(check.dist(list(x=Inf,fx=1)))
})
