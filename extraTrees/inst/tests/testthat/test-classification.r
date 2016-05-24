
context("Classification basics")

get.data.class <- function(n=800) {
  p = 4
  x = matrix(runif(n*p), n, p)
  f = function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[,3]>0.4)
  y = factor(f(x), levels=0:3)
  return(list(x=x, y=y))
}

train <- get.data.class(100)
test  <- get.data.class(101)

test_that("basic classification and prediction", {  
  et   <- extraTrees(train$x, train$y, nodesize=1, mtry=2, numRandomCuts=2, ntree=50)
  expect_equal( 50, et$ntree )
  expect_true(  et$factor )  
  expect_false( et$multitask )

  ## standard prediction
  yhat <- predict(et, test$x)
  expect_true ( is.factor(yhat) )
  expect_equal( length(yhat), length(test$y) )
  expect_equal( levels(yhat), levels(test$y) )
  
  ## allValues prediction
  yall = predict(et, test$x, allValues=T)
  expect_true ( is.factor(yall[,1]) )
  expect_equal( nrow(yall), nrow(test$x) )
  expect_equal( ncol(yall), 50 )
  expect_false( is.numeric(yall) )
  
  ## probability prediction
  yprob = predict(et, test$x, probability=T)
  expect_true ( is.numeric(yprob) )
  expect_equal( nrow(yprob), nrow(test$x) )
  expect_equal( ncol(yprob), nlevels(train$y) )
  expect_equal( colnames(yprob), levels(train$y) )  
  expect_equal( unname(mean(yall[1,]==0)), unname(yprob[1,1]), tolerance=1e-6 )
  expect_equal( rowSums(yprob), rep.int(1, nrow(yprob)) )
})

test_that("selectTree works with classification", {  
  et    <- extraTrees(train$x, train$y, nodesize=1, mtry=2, numRandomCuts=2, ntree=20)
  trees <- rep(c(T, F), 10)
  et10 <- selectTrees(et, trees )
  expect_equal( et10$ntree, 10 )
  yall   <- predict(et,   test$x, allValues=T)
  yall10 <- predict(et10, test$x, allValues=T)
  expect_true( all(yall[,trees] == yall10) )
})
