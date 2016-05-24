context("kFoldCV")

test_that("various k and .chunkSize can be inputted", {
  f <- function(data, newdata, param1, param2) {
    if ( (param1 * param2) %% 2 == 0) {
      param1 + 0.01
    } else {
      param2 - 0.01
    }
  }
  
  parameterGrid <- expand.grid(param1=1:3, param2=1:3)
  
  # function is independent of train and test
  correctOutput <- foreach(row=iter(parameterGrid, by="row"), .combine=c) %do% {
    g <- function(...) f(data=NULL, new=NULL, ...)
    do.call(g, as.list(row))
  }
  
  # mix and match k and .chunkSize parameters
  testingGrid <- expand.grid(k=c(3,5,7),.chunkSize=c(1,5,10))
  
  foreach(testCase=iter(testingGrid, by="row")) %do% {
    g <- function(...) kFoldCV(proc=f, data=USArrests, params=parameterGrid, ...)
    
    suppressWarnings(testOutput <- do.call(g, as.list(testCase)))
    expect_that(testOutput, equals(correctOutput))
  }
  
  
})
test_that("results are reproducible", {
  require(class)
  f <- function(data, newdata, k) {
    preds <- knn(train=data[,-1], test=newdata[, -1], cl=data[, 1], k=k)
    mean(preds==newdata[, 1])
  }
  .iris <- iris[, 5:1]
  parameterGrid <- expand.grid(k=3:5)
  out1 <- kFoldCV(f, 10, .iris, parameterGrid, .rngSeed=407)
  out2 <- kFoldCV(f, 10, .iris, parameterGrid, .rngSeed=407)
  expect_that(out1, equals(out2))
})

test_that("kFoldCV doc examples work", {
    # simple example with k-NN where we can build our own wrapper
    library(class)
    data(iris)
    .iris <- iris[, 5:1] # put response as first column
    
    # make a wrapper for class::knn
    f <- function(data, newdata, k) {
      preds <- knn(train=data[,-1],
                   test=newdata[, -1],
                   cl=data[, 1],
                   k=k)
      mean(preds==newdata[, 1])
    }
    
    params <- list(k=seq.int(3))
    
    accuracy <- kFoldCV(f, 10, .iris, params, .rngSeed=407)
    
    expect_that(is.null(data.frame(expand.grid(params), accuracy=accuracy)),
                is_false())
    
    # look at a more complicated example:
    # cross validate an svm with different kernels and different models
    require(e1071)
    g <- function(data, newdata, kernel, cost, gamma, formula) {
      kern <- switch(kernel, "linear", "radial", stop("invalid kernel"))
      form <- switch(formula,
                     as.formula(Species ~ .),
                     as.formula(Species ~ Petal.Length + Petal.Width),
                     as.formula(Petal.Length ~ .),
                     stop('invalid formula'))
      
      svmWrapper <- function(data, newdata, kernel, cost, gamma, form) {
        svmObj <- svm(formula=form, data=data, kernel=kernel,
                      cost=cost, gamma=gamma)
        predict(svmObj, newdata)
      }
      
      preds <- svmWrapper(data, newdata, kernel=kern, cost=cost,
                          gamma=gamma, form=form)
      
      if (formula != 3) {
        mean(preds == newdata[["Species"]])
      } else {
        mean((preds - newdata[["Petal.Length"]])^2)
      }
    }
    
    params <- list(kernel=1:2, cost=c(10,50), gamma=0.01, formula=1)
    accuracy <- kFoldCV(g, 10, iris, params)
    
    expect_that(is.null(data.frame(expand.grid(params), metric=accuracy)),
                is_false())
    
})