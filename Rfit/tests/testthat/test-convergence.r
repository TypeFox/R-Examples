EPS <- sqrt(.Machine$double.eps)

test_that("rfit convergence", {
  x<-y<-rep(0:1,each=10)
  err<-rfit(y~x)$coef[2]-1
  expect_less_than(err,EPS)
})


