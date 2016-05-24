
context("Quantile regression")

test_that("quantile regression and prediction", {
  n = 400
  ntest = 200
  p = 4
  x = matrix(runif(n*p), n, p)
  f = function(x) (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
  y = as.numeric(f(x))
  
  xtest = matrix(runif(ntest*p), ntest, p)
  ytest = f(xtest)
  
  ## learning extra trees:
  et = extraTrees(x, y, mtry=p, quantile=T, ntree=100)
  expect_is   ( et, "extraTrees")
  expect_false( et$factor )
  expect_false( et$multitask )
  expect_true(  et$quantile  )
  expect_equal( 100, et$ntree )
  
  yhat0.5 = predict(et, xtest, quantile = 0.5)
  yhat0.8 = predict(et, xtest, quantile = 0.8)
  expect_equal(length(yhat0.5), nrow(xtest))
  expect_true(all(yhat0.5 <= yhat0.8))
})
