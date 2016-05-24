context("constants for other unit tests")

suppressPackageStartupMessages(library(class))
suppressPackageStartupMessages(library(e1071))

df <- within(iris[40:60, ], {
  Setosa <- as.factor(2*as.numeric(Species == "setosa")-1)
  Species <- NULL
})

preds <- iris[sample(nrow(iris), size=15), 5]

response <- iris[seq.int(15), 5]

wts <- runif(n=15, min=0, max=1)


#########################################
# contrived test procedures for unit tests
#########################################

testSVM <- buildEstimationProcedure(train=svm)

testSVMProcArgs <- list(.trainArgs=list(formula= Setosa ~ ., cost=10))


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

##################################
# contrived reweighter for unit tests.
##################################

testReweighterMetadata <- list(reweighterInputWts="w", reweighterInputResponse="truth",
                            reweighterInputPreds="preds", reweighterOutputWts="w")

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

######################################
# contrived aggregator for unit tests
######################################

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

##################################
# contrived performance analyzers
##################################
testPerfAnalyzer <- function(prediction, response, oobObs, Xx) { 
  c(x=Xx)
}

testPerfAnalyzer2 <- function(pred, truth, oob, zeta) {
  list(e=mean(pred != truth), z=zeta)
}

testPerfAnalyzer2Metadata <- list(analyzerInputPreds="pred",
                                  analyzerInputResponse="truth",
                                  analyzerInputOOBObs = "oob")