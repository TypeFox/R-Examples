require(LiblineaR)
require(SwarmSVM)

context("clusterSVM")

data(svmguide1)
svmguide1.t = svmguide1[[2]]
svmguide1 = svmguide1[[1]]

data(iris)

test_that("Error Trigger",{
  # Wrong parameters
  expect_error({csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 0,
                                     centers = 8, seed = 512, verbose = 0,
                                     valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  })
  expect_error({csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                                      centers = 8, seed = 512, verbose = 0, type = 11,
                                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  })
  expect_error({csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                                      centers = 8, seed = 512, verbose = 0, cost = -1,
                                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  })
  expect_error({csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                                      centers = 8, seed = 512, verbose = 0, epsilon = -1,
                                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  })
  expect_error({csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                                      centers = 8, seed = 512, verbose = 0, verbose = 10,
                                      valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  })
  
  # Wrong prediction object
  csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                        centers = 8, seed = 512, verbose = 0,
                        valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  expect_error({pred = predict(csvm.obj$sparse, svmguide1.t[,-1])})
  
  # Only one cluster
  expect_warning({
    csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                          centers = 1, seed = 512, verbose = 0,
                          valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  })
})

test_that("Switch Clustering function",{
  csvm.obj.1 = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                          centers = 8, seed = 512, verbose = 0, 
                          valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1],
                          cluster.method = "kmeans")
  csvm.obj.2 = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                          centers = 8, seed = 512, verbose = 0, 
                          valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1],
                          cluster.method = "mlKmeans")
  # Avoid the error from kkmeans
  csvm.obj.3 = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                          centers = 2, seed = 1024, verbose = 0, 
                          valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1],
                          cluster.method = "kernkmeans")
  expect_true(csvm.obj.3$time$total.time>csvm.obj.2$time$total.time)
  expect_true(csvm.obj.3$time$total.time>csvm.obj.1$time$total.time)
})

test_that("Performance",{
  liblinear.obj = LiblineaR::LiblineaR(data = svmguide1[,-1], target = svmguide1[,1], 
                                       type = 1, verbose = FALSE)
  liblinear.pred = predict(liblinear.obj, svmguide1.t[,-1])$prediction
  liblinear.score = sum(liblinear.pred==svmguide1.t[,1])/length(liblinear.pred)
  
  # with 8 clustering
  csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                        centers = 8, seed = 512, verbose = 0,
                        valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  csvm.score = csvm.obj$valid.score
  expect_true(csvm.score>liblinear.score)
  
  # with 1 clustering
  expect_warning({
    csvm.obj = clusterSVM(x = svmguide1[,-1], y = svmguide1[,1], lambda = 1,
                          centers = 1, seed = 512, verbose = 0,
                          valid.x = svmguide1.t[,-1],valid.y = svmguide1.t[,1])
  })
  csvm.score = csvm.obj$valid.score
  expect_equal(csvm.score, liblinear.score)
  
  # Multiclassification
  csvm.obj = clusterSVM(x = as.matrix(iris[,-5]), y = iris[,5], sparse = FALSE,
                        centers = 2, seed = 512, verbose = 0,
                        valid.x = as.matrix(iris[,-5]),valid.y = iris[,5])
  expect_true(csvm.obj$valid.score>0.97)
})

set.seed(512)
xorx = rbind(cbind(runif(100,1,2),runif(100,1,2)),
             cbind(runif(100,-2,-1),runif(100,1,2)),
             cbind(runif(100,1,2),runif(100,-2,-1)),
             cbind(runif(100,-2,-1),runif(100,-2,-1)))
xory = c(rep(1,100),
         rep(0,200),
         rep(1,100))


test_that("XOR Toy Data",{
  
  liblinear.obj = LiblineaR::LiblineaR(data = xorx, target = xory, 
                                       type = 1, verbose = FALSE)
  liblinear.pred = predict(liblinear.obj, xorx)$prediction
  liblinear.score = sum(liblinear.pred==xory)/length(liblinear.pred)
  expect_true(liblinear.score<0.6)
  
  csvm.obj = clusterSVM(x = xorx, y = xory, lambda = 1,
                        centers = 2, seed = 512, verbose = 0,
                        valid.x = xorx, valid.y = xory,
                        cluster.method = 'kmeans')
  expect_equal(csvm.obj$valid.score,1)
  
  csvm.obj = clusterSVM(x = xorx, y = xory, lambda = 1,
                        centers = 2, seed = 512, verbose = 0,
                        valid.x = xorx, valid.y = xory,
                        cluster.method = 'mlKmeans')
  expect_equal(csvm.obj$valid.score,1)
  
  csvm.obj = clusterSVM(x = xorx, y = xory, lambda = 1,
                        centers = 2, seed = 1024, verbose = 0,
                        valid.x = xorx, valid.y = xory,
                        cluster.method = 'kernkmeans')
  expect_equal(csvm.obj$valid.score,1)
})
