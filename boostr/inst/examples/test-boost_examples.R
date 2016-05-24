### Demonstrate simple call with just list(train=svm)

library(foreach)
library(iterators)
library(e1071)

svmArgs <- list(formula=Species~., cost=100)
boost(x=list(train=svm),
      reweighter=arcfsReweighter,
      aggregator=arcfsAggregator,
      data=iris,
      .procArgs=list(.trainArgs=svmArgs),
      B=2)

### Demonstrate call with train and predict and custom 
### reweighters and aggregators

df <- within(iris, {
  Setosa <- as.factor(2*as.numeric(Species == "setosa")-1)
  Species <- NULL
})

# custom predict function
newPred <- function(obj, new) {
  predict(obj, new)
}

predMetadata <- c(modelName="obj", predictionSet="new")

# custom reweighter
testReweighterMetadata <- list(
                            reweighterInputWts="w",
                            reweighterInputResponse="truth",
                            reweighterInputPreds="preds",
                            reweighterOutputWts="w")

testReweighter <- function(preds, truth, w) {
  
  wrongPreds <- (preds != truth)
  err <- mean(wrongPreds)
  if (err != 0) {
    new_w <- w / err^(!wrongPreds)
  } else {
    new_w <- runif(n=length(w), min=0, max=1)
  }
  
  
  list(w=new_w, alpha=rnorm(1))
}

# custom aggregator
testAggregatorMetadata <- c(.inputEnsemble="ensemble")

testAggregator <- function(ensemble) {
  weights <- runif(min=0, max=1, n=length(ensemble))
  function(x) {
    preds <- foreach(estimator = iter(ensemble),
                     .combine = rbind) %do% {
                       matrix(as.character(estimator(x)), nrow=1)
                     }
    
    as.factor(predictClassFromWeightedVote(preds, weights))
  }
}

# collect all the relevant metadata
metadata <- c(predMetadata, testReweighterMetadata, testAggregatorMetadata)

# set additional procedure arguments
procArgs <- list(
              .trainArgs=list(
                formula=Setosa ~ .,
                cost=100)
              )

#test boost when irrelevant metadata is passed in.
boostedSVM <- boost(list(train=svm, predict=newPred),
                    B=3,
                    reweighter=testReweighter,
                    aggregator=testAggregator,
                    data=df,
                    metadata=metadata,
                    .procArgs=procArgs,
                    .boostBackendArgs=list(
                      .reweighterArgs=list(fakeStuff=77))
                    )

### Demonstrate customizing 'metadata' for estimation procedure
library(class)

testkNNProcMetadata <- list(learningSet="traindata", predictionSet="testdata")

testkNNProc <- function(formula, traindata, k) {  
  df <- model.frame(formula=formula, data=traindata)
  function(testdata, prob=FALSE) {
    df2 <- tryCatch(model.frame(formula=formula, data=testdata)[, -1],
                    error = function(e) testdata 
    )
    knn(train=df[, -1], test=df2, cl=df[, 1], prob=prob, k=k) 
  }
}

testKNNProcArgs <- list(formula=Setosa ~ ., k = 5)

metadata <- testkNNProcMetadata
boostBackendArgs <- list(.reweighterArgs=list(m=0))

boostedKNN <- boost(x=testkNNProc, B=3,
      reweighter=arcx4Reweighter,
      aggregator=arcx4Aggregator,
      data=df, 
      metadata=metadata,
      .boostBackendArgs=boostBackendArgs,
      .procArgs=testKNNProcArgs)

### Demonstrate using an alternative performance analyzer

testPerfAnalyzer2 <- function(pred, truth, oob, zeta) {
  list(e=mean(pred != truth), z=zeta)
}

testPerfAnalyzer2Metadata <- list(analyzerInputPreds="pred",
                                  analyzerInputResponse="truth",
                                  analyzerInputOObObs="oob")

metadata <- c(metadata, testPerfAnalyzer2Metadata)

boostedkNN <- boost(testkNNProc,
                    B=3,
                    reweighter=vanillaBagger,
                    aggregator=vanillaAggregator,
                    data=df,
                    .procArgs=testKNNProcArgs,
                    metadata=metadata,
                    .boostBackendArgs = list(
                      .analyzePerformanceArgs = list(zeta="77"),
                      .reweighterArgs=list(fakeStuff=77)),
                    analyzePerformance=testPerfAnalyzer2)

