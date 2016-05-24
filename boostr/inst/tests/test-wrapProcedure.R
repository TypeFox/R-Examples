context("wrapProcedure")

library(class)

# an estimation procedure outputs a estimator

knnProc <- function(formula, traindata, k) {
  df <- model.frame(formula=formula, data=traindata)
  function(testdata, prob=FALSE) {
    knn(train=df[, -1], test=testdata, cl=df[, 1], prob=prob, k=k) 
  }
}

boostrKNN <- wrapProcedure(knnProc,
                          learningSet="traindata",
                          predictionSet="testdata")

df <- model.frame(formula=Species ~ ., data=iris)[1:120, ]

boostrEstimator <- do.call(boostrKNN, list(data=df, k=5, formula=Species~.))

test_that("wrapProcedure's output adheres to spec. in whitepaper", {
  expect_that(boostrKNN, is_a('function'))
  
  factoryArgs <- names(formals(boostrKNN))
  expect_that(mean(c("data") %in% factoryArgs), equals(1))
  

  expect_that(boostrEstimator, is_a('function'))
  
  estimatorArgs <- names(formals(boostrEstimator))
  expect_that(c("newdata") %in% estimatorArgs, is_true())
})

test_that("wrapProcedure's output doesn't affect estimation", {
  expect_that(boostrEstimator(iris[120:130,-5], .estimatorArgs=list(prob=TRUE)),
              equals(knnProc(Species ~ . , df, k=5)(iris[120:130, -5], TRUE))
  )
})

test_that("wrapProcedure's output is compatible with boostBackend", {
  phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(iris)),
                reweighter=arcx4Reweighter,
                aggregator=arcx4Aggregator,
                proc=boostrKNN,
                .procArgs=list(formula=Species~., k=5),
                .reweighterArgs=list(m=0),
                data=iris,
                .formatData=TRUE)
  
  expect_that(is.null(phi), is_false())
})

