library("fastAdaboost")
test_that("rpart call from Rcpp works",{
  num_each <- 100
  fakedata <- data.frame( X=c(rnorm(num_each,0,1),rnorm(num_each,1.5,1)), Y=c(rep(0,num_each),rep(1,num_each) ) )
  fakedata$Y <- factor(fakedata$Y)
  num_examples <- nrow(fakedata)
  weights <- rep(1./num_examples, num_examples)
  classname_map <- levels(fakedata[,"Y"])
  names(classname_map) = c("0","1")
  x <- call_rpart_(Y~X,wrap_rpart, fakedata, weights, classname_map)
  
  tree <- x$tree
  predictions <- predict(tree, newdata=fakedata, type="class")
  #print(predictions)
  
})