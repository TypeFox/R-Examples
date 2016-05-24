
context("Weights")

test_that("weighted regression works", {
  n = 200
  ntest = 100
  p = 4
  x = matrix(runif(n*p), n, p)
  f = function(x) (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
  y = as.numeric(f(x))
  weights = rep.int(c(1, 0.5), n/2)
  
  xtest = matrix(runif(ntest*p), ntest, p)
  ytest = f(xtest)
  
  ## learning extra trees:
  et = extraTrees(x, y, mtry=p, ntree=100, weights=weights)
  expect_is   ( et, "extraTrees")
  expect_false( et$factor )
  expect_false( et$multitask )
  expect_true( et$useWeights )
  
  yhat = predict(et, xtest)
  expect_equal(length(yhat), nrow(xtest))
})
