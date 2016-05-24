library("fastAdaboost")
test_that("running adaboost works",{
  #create data
  num_each <- 1000
  fakedata <- data.frame( X=c(rnorm(num_each,0,1),rnorm(num_each,1.5,1)), Y=c(rep(0,num_each),rep(1,num_each) ) )
  fakedata$Y <- factor(fakedata$Y)
  #run adaboost
  A <- adaboost(Y~X, fakedata, 10)
})