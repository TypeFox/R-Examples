
context("Multi-task learning")

get.data.class <- function(n=800) {
  p = 4
  x = matrix(runif(n*p), n, p)
  task = sample(1:2, size=n, replace=T)
  f = function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[cbind(1:n,task+2)]>0.4)
  y = factor(f(x) %% 2, levels=0:1)
  yd = f(x) + 0.2*rnorm(n)
  return(list(x=x, y=y, yd=yd, task=task))
}

train <- get.data.class(100)
test  <- get.data.class(101)

test_that("mt classification and prediction", {  
  et   <- extraTrees(train$x, train$y, tasks=train$task, numRandomCuts=2, ntree=50)
  expect_equal( 50, et$ntree )
  expect_true(  et$factor )  
  expect_true( et$multitask )
  
  ## prediction
  expect_error( predict(et, test$x) )
  
  yhat <- predict(et, test$x, newtasks=test$task)
  expect_equal( length(yhat), length(test$y) )
  expect_equal( levels(yhat), levels(test$y) )
  
  ## allValues prediction
  expect_error( predict(et, test$x, allValues=T) )
  yall = predict(et, test$x, newtasks=test$task, allValues=T)
  expect_equal( nrow(yall), nrow(test$x) )
  expect_equal( ncol(yall), 50 )
  expect_false( is.numeric(yall) )
  expect_true ( is.factor(yall[,1]) )
})

test_that("mt regression and prediction", {  
  et   <- extraTrees(train$x, train$yd, tasks=train$task, numRandomCuts=2, ntree=50)
  expect_equal( 50, et$ntree )
  expect_false( et$factor   )
  expect_true( et$multitask )
  
  ## prediction
  expect_error( predict(et, test$x) )
  
  yhat <- predict(et, test$x, newtasks=test$task)
  
  expect_equal( length(yhat), length(test$yd) )
  
  ## allValues prediction
  expect_error( predict(et, test$x, allValues=T) )
  yall = predict(et, test$x, newtasks=test$task, allValues=T)
  expect_equal( nrow(yall), nrow(test$x) )
  expect_equal( ncol(yall), 50 )
  expect_true( is.numeric(yall) )
  expect_equal( yhat, rowMeans(yall), tolerance=1e-5)
})

test_that("mt classification with more than 2 classes fails", {
  y = factor(rep_len(1:5, nrow(train$x)))
  expect_error( extraTrees(train$x, y, tasks = train$task), 
                regexp = "Multi-task learning only works with 2 factors" )
})
