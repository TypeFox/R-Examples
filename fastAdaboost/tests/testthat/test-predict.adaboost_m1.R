library("fastAdaboost")
test_that("predicting adaboost works",{
  #create data
  num_each <- 1000
  fakedata <- data.frame( X=c(rnorm(num_each,0,1),rnorm(num_each,1.5,1)), Y=c(rep(0,num_each),rep(1,num_each) ) )
  fakedata$Y <- factor(fakedata$Y)
  #run adaboost
  A <- adaboost(Y~X, fakedata, 10)
  #print(A)
  pred <- predict(A,newdata=fakedata)
  print(paste("Adaboost Error on fakedata:",pred$error))
  print(table( pred$class, fakedata$Y))
  expect_true(pred$error<1.)
})

test_that("adaboost errors are predicted correctly",{
  num_each <- 1000
  fakedata <- data.frame( X=c(rnorm(num_each,0,1),rnorm(num_each,1.5,1)), Y=c(rep(0,num_each),rep(1,num_each) ) )
  fakedata$Y <- factor(fakedata$Y)
  #run adaboost
  A <- adaboost(Y~X, fakedata, 10)
  pred <- predict(A,newdata=fakedata)
  err <- length(which(pred$class!=fakedata$Y))/nrow(fakedata)
  expect_true( abs(err - pred$error)<1e-5 ) 
})
