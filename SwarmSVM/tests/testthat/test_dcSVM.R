require(SwarmSVM)

context("dcSVM")

data(svmguide1)
svmguide1.t = svmguide1[[2]]
svmguide1 = svmguide1[[1]]

data(iris)

test_that("Error Trigger",{
  expect_error({
    model = dcSVM(x = as.matrix(iris[,-5]), y = iris[, 5], m = 100, k = 2, max.levels = 1)
  })
  
  expect_error({
    model = dcSVM(x = as.matrix(iris[,-5]), y = iris[, 5], m = 100, max.levels = 1, k = 100)
  })
  
  expect_error({
    model = dcSVM(x = as.matrix(iris[,-5]), y = iris[, 5], m = 100, k = 1, max.levels = 1)
  })
  
})

test_that("Performance",{
  dcsvm.model = dcSVM(x = svmguide1[,-1], y = svmguide1[,1],
                      k = 4, max.levels = 4, seed = 0, cost = 32, gamma = 2,
                      kernel = 3,early = 0, m = 800, pre.scale = TRUE,
                      valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1])
  score = dcsvm.model$valid.score
  expect_true({score>0.9})
  
  dcsvm.model = dcSVM(x = as.matrix(svmguide1[,-1]), y = svmguide1[,1], k = 10, max.levels = 1, 
                      early = 1, gamma = 2, cost = 32, tolerance = 1e-2, m = 800, 
                      pre.scale = TRUE,
                      valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1])
  score = dcsvm.model$valid.score
  expect_true({score>0.9})
})

test_that("Sparse Input",{
  local.file.name = tempfile()
  download.file("http://www.sfu.ca/~hetongh/data/ijcnn1.dcsvm.RData",local.file.name)
  load(local.file.name)
  ijcnn1.t = ijcnn1[[2]]
  ijcnn1 = ijcnn1[[1]]
  
  dcsvm.model = dcSVM(x = ijcnn1[,-1], y = ijcnn1[,1], k = 10, max.levels = 1, seed = 1024,
                      early = 1, gamma = 2, cost = 32, tolerance = 1e-2, m = 5000, scale = FALSE,
                      valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
  
  expect_error({dcsvm.model = dcSVM(x = as.matrix(ijcnn1[,-1]), y = ijcnn1[,1], k = 10, max.levels = 1, seed = 1024,
                                    early = 1, gamma = 2, cost = 32, tolerance = 1e-2, m = 5000, scale = FALSE,
                                    valid.x = ijcnn1.t[,-1], valid.y = ijcnn1.t[,1])
  })
  
})
