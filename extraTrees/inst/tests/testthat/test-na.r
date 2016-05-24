context("Missing x values")

get.data.class <- function(n=800) {
  p = 4
  x = matrix(runif(n*p), n, p)
  f = function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[,3]>0.4)
  y = factor(f(x), levels=0:3)
  x[1, 1] = NA
  x[2, 2] = NA
  x[3, 3] = NA
  x[4, 4] = NA
  return(list(x=x, y=y))
}

train <- get.data.class(100)
test  <- get.data.class(101)

test_that("na.action='stop' gives error", {  
  expect_error( extraTrees(train$x, train$y, ntree=50, na.action='stop'), regex="NA" )
})

test_that("na.action='zero' does not give error", {  
  et   <- extraTrees(train$x, train$y, ntree=50, na.action='zero')
  yhat <- predict(et, test$x)
  expect_equal(et$xHasNA, FALSE)
})

test_that("na.action='fuse' enables fusing", {  
  et   <- extraTrees(train$x, train$y, ntree=50, na.action='fuse')
  yhat <- predict(et, test$x)
  expect_equal(et$xHasNA, TRUE)
})
