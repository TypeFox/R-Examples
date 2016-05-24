
library("Matrix")

context("Sparse classification")

test_that("sparse matrix gets all values", {
  s1 <- sparseMatrix( i=1:4, j=2:5, x=c(1,-2,0,-1) )
  js1 <- toJavaCSMatrix( s1 )
  expect_equal( .jcall(js1, "I", "getNumNonZero"), 3 )
})


get.sparse.class <- function(n=800) {
  p = 40
  x = matrix(sample(0:1, n*p, prob=c(0.9, 0.1), replace=T), n, p)
  x = as(x, "sparseMatrix")
  f = function(x) (Matrix::rowSums(x[,1:20]) >= 3) + (Matrix::rowSums(x[,21:p]) >= 4)
  y = factor(f(x), levels=0:2)
  return(list(x=x, y=y))
}
 
train <- get.sparse.class(100)
test  <- get.sparse.class(101)

test_that("classification using sparse input", {  
  et   <- extraTrees(train$x, train$y, nodesize=1, ntree=50)
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
})
