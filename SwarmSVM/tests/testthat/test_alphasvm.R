require(SwarmSVM)

context("alphasvm")

data(svmguide1)
svmguide1.t = svmguide1[[2]]
svmguide1 = svmguide1[[1]]

data(iris)

test_that("Error Trigger",{
  expect_error({ model = alphasvm(x = iris[,-5], y = iris[,5], alpha = 0) })
  expect_error({ model = alphasvm(x = iris[,-5], y = iris[,5], alpha = rep(0,nrow(iris))) })
  expect_error({ model = alphasvm(x = svmguide1[,-1], y = svmguide1[1:100,1])})
  model = alphasvm(x = iris[,-5], y = iris[,5])
  expect_error({preds = predict(model, iris[,1:3])})
  expect_error({preds = predict(1, iris[,-5])})
  preds = predict(model, iris[,-5])
})

test_that("Performance",{
  set.seed(1024)
  time.stamp = proc.time()
  model = alphasvm(x = svmguide1[,-1], y = svmguide1[,1], cost = 32)
  time.elapse1 = (proc.time()-time.stamp)[3]
  preds = predict(model, svmguide1.t[,-1])
  score = sum(diag(table(preds,svmguide1.t[,1])))/nrow(svmguide1.t)
  expect_true(score>0.8)
  
  # Take the previous alpha
  new.alpha = matrix(0, nrow(svmguide1),1)
  new.alpha[model$index,] = model$coefs
  time.stamp = proc.time()
  model2 = alphasvm(x = svmguide1[,-1], y = svmguide1[,1], alpha = new.alpha, cost = 32)
  time.elapse2 = (proc.time()-time.stamp)[3]
  preds = predict(model2, svmguide1.t[,-1])
  score = sum(diag(table(preds,svmguide1.t[,1])))/nrow(svmguide1.t)
  expect_true(score>0.8)
  
  # Take random alpha
  new.alpha = matrix(rnorm(nrow(svmguide1)), nrow(svmguide1),1)
  new.alpha[model$index,] = model$coefs
  time.stamp = proc.time()
  model2 = alphasvm(x = svmguide1[,-1], y = svmguide1[,1], alpha = new.alpha, cost = 32)
  time.elapse3 = (proc.time()-time.stamp)[3]
  preds = predict(model2, svmguide1.t[,-1])
  score = sum(diag(table(preds,svmguide1.t[,1])))/nrow(svmguide1.t)
  expect_true(score>0.8)
  
  # Running time is not so consistent
  # expect_true(time.elapse3 > time.elapse1 && time.elapse1 > time.elapse2)
  expect_true(time.elapse3 > time.elapse2 && time.elapse1 > time.elapse2)
})


