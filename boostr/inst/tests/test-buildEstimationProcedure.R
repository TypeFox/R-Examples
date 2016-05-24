context("buildEstimationProcedure")

test_that("buildEstimationProcedure agrees with hand-made 'lm' wrapper", {
  form <- formula(Petal.Length ~ Petal.Width)
  
  Phi <- buildEstimationProcedure(lm)
  
  set.seed(1234)
  ibIndx <- sample(nrow(iris), prob=rep.int(1, nrow(iris)), replace=TRUE)
  resampledData <- iris[ibIndx, ]
  
  phi <- Phi(data=resampledData, .trainArgs = list(formula=form))
  
  
  lmObj <- lm(formula=form, data=resampledData)
  phi2 <- function(newdata) predict(lmObj, newdata)
  
  expect_that(phi(NULL), equals(phi2(NULL)))
})

suppressPackageStartupMessages(require(randomForest))

Phi <- buildEstimationProcedure(train=randomForest)
phi <- Phi(data=iris, .trainArgs=list(formula=Species ~ .))

test_that("buildEstimationProcedure yields boostr compatible wrappers", {

  # check the estimation procedure has the proper signature
  test1 <- mean(names(formals(Phi)) %in%
              c("data", "weights", ".trainArgs", ".predictArgs"))
  expect_that(test1, equals(1))
  
  # check the estimator has the proper signature
  test2 <- mean(names(formals(phi)) %in% c("newdata", ".predictArgs"))
  expect_that(test2, equals(1))
})

test_that("resulting estimator doesn't goof and can accept .predictArgs", {
  preds <- phi(iris[1:15, ])
  
  expect_that(length(preds), equals(15))
  expect_that(preds, is_a("factor"))

  probs <- phi(iris[1:15, ], list(type='prob'))
  
  expect_that(dim(probs), equals(c(15,3)))
  expect_that(rowSums(probs), is_equivalent_to(rep.int(1, nrow(probs))))
})

test_that("buildEstimationProcedure is compatible with boostBackend", {
  procArgs <- list(.trainArgs=list(formula=Species ~ . ))
  phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(iris)),
                reweighter=arcx4Reweighter,
                aggregator=arcx4Aggregator,
                proc=Phi,
                .procArgs=procArgs,
                data=iris,
                .reweighterArgs=list(m=0),
                .formatData=TRUE)
  
  expect_that(is.null(phi), is_false())
})