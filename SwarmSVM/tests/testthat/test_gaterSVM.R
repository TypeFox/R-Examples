require(SwarmSVM)

context("gaterSVM")

data(svmguide1)
svmguide1.t = svmguide1[[2]]
svmguide1 = svmguide1[[1]]

# Toy example 
set.seed(1024)
train.1 = cbind(runif(333,-1.7,-0.7),
                runif(333,0.7,1.7))
train.2 = cbind(runif(333,-0.5,0.5),
                runif(333,-0.5,0.5))
train.3 = cbind(runif(334,0.7,1.7),
                runif(334,-1.7,-0.7))
y = c(rep(1,333),rep(0,333),rep(1,334))
train = rbind(train.1, train.2, train.3)
train = cbind(train, y)

test.1 = cbind(runif(3333,-1.7,-0.7),
               runif(3333,0.7,1.7))
test.2 = cbind(runif(3333,-0.5,0.5),
               runif(3333,-0.5,0.5))
test.3 = cbind(runif(3334,0.7,1.7),
               runif(3334,-1.7,-0.7))
y = c(rep(1,3333), rep(0,3333), rep(1,3334))
test = rbind(test.1, test.2, test.3)
test = cbind(test, y)

toydata = list(train, test)
toydata.t = test
toydata = train

test_that("Error Trigger",{
  expect_error({
    gaterSVM.model = gaterSVM(x = svmguide1[,-1], y = svmguide1[,1], hidden = -1, seed = 0,
                              m = 10, max.iter = 3, learningrate = 0.01, threshold = 1, stepmax = 1000,
                              valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1], verbose = TRUE)
  })
  expect_error({
    gaterSVM.model = gaterSVM(x = svmguide1[,-1], y = svmguide1[,1], hidden = 10, seed = 0,
                              m = 10, max.iter = 0, learningrate = 0.01, threshold = 1, stepmax = 1000,
                              valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1], verbose = TRUE)
  })
  expect_error({
    gaterSVM.model = gaterSVM(x = svmguide1[,-1], y = svmguide1[,1], hidden = 10, seed = 0,
                              m = 0, max.iter = 3, learningrate = 0.01, threshold = 1, stepmax = 1000,
                              valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1], verbose = TRUE)
  })
  expect_error({
    gaterSVM.model = gaterSVM(x = svmguide1[,-1], y = svmguide1[,1], hidden = 10, seed = 0, c = -1,
                              m = 10, max.iter = 3, learningrate = 0.01, threshold = 1, stepmax = 1000,
                              valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1], verbose = TRUE)
  })
})

test_that("Performance",{
  gaterSVM.model = gaterSVM(x = toydata[,-3], y = toydata[,3], hidden = 10, seed = 0,
                            m = 5, max.iter = 1, learningrate = 0.001, threshold = 0.05,
                            valid.x = toydata.t[,-3], valid.y = toydata.t[,3], verbose = TRUE)
  expect_true(gaterSVM.model$valid.score == 1)
  gaterSVM.model = gaterSVM(x = svmguide1[,-1], y = svmguide1[,1], hidden = 10, seed = 0,
                            m = 10, max.iter = 3, learningrate = 0.01, threshold = 1, stepmax = 1000,
                            valid.x = svmguide1.t[,-1], valid.y = svmguide1.t[,1], verbose = TRUE)
  expect_true(gaterSVM.model$valid.score>0.9)
})
