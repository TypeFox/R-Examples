
## ----, kNNExample, cache=FALSE-------------------------------------------
kNN_EstimationProcedure <- function(k, learningSet) {
  learningInput <- learningSet[, -1]
  learningResponse <- learningSet[, 1]
  
  function(newdata) {
    class::knn(train = learningInput,
               cl = learningResponse,
               k = k,
               test = newdata)  
  }
}

library(mlbench)
data(Glass)

# train estimator on Glass dataset
kNN_Estimator <- kNN_EstimationProcedure(5, Glass[,10:1])

# predict first 10 observations of Glass data
kNN_Estimator(Glass[1:10, 9:1])

table(kNN_Estimator(Glass[1:10, 9:1]), Glass[1:10,10])


## ----, svmExample, cache=FALSE-------------------------------------------
svm_EstimationProcedure <- function(formula, cost, data) {
  model <- e1071::svm(formula, cost=cost, data=data)  
  function(newdata) {
    predict(model, newdata=newdata)
  }
}

# train estimator on Glass dataset
svm_Estimator <- svm_EstimationProcedure(formula(Type ~ .), 100, Glass)

# predict first 10 observations of Glass data
svm_Estimator(Glass[1:10,-10])

table(svm_Estimator(Glass[1:10,-10]), Glass[1:10,10])


## ----, helpFxnTable, echo=FALSE, cache=FALSE, results='asis'-------------
helperFxns <- matrix(c("wrapProcedure", "buildEstimationProcedure", "wrapReweighter", "wrapAggregator", "WrapPerformanceAnalyzer"), ncol=1)
colnames(helperFxns) <- "Wrapper Generators"

helpFxnXTab <- xtable::xtable(helperFxns, caption="Table of `boostr`'s Wrapper Generators")

print(helpFxnXTab, include.rownames=FALSE, type="html")


## ----, arcx4AndkNN, cache=FALSE------------------------------------------
boostr::boostWithArcX4(x = kNN_EstimationProcedure,
                       B = 3,
                       data = Glass,
                       metadata = list(learningSet="learningSet"),
                       .procArgs = list(k=5),
                       .boostBackendArgs = list(
                         .subsetFormula=formula(Type~.))
                       ) 


## ----, arcx4AndSvm, cache=FALSE------------------------------------------
boostr::boostWithArcX4(x = list(train = e1071::svm),
                       B = 3,
                       data = Glass,
                       .procArgs = list(
                         .trainArgs=list(
                           formula=formula(Type~.),
                           cost=100
                           )
                         )
                       )


## ----, cache=FALSE-------------------------------------------------------
boostr::arcx4Reweighter


## ----, reweighterExample, cache=FALSE------------------------------------
exoticReweighter <- function(wts, truth, preds) {
  permutedWts <- sample(wts)
  list(wts=permutedWts)
}

boostr::boost(x = list(train=e1071::svm), B = 3,
              initialWeights = seq.int(nrow(Glass)),
              reweighter = exoticReweighter,
              aggregator = boostr::vanillaAggregator,
              data = Glass,
              .procArgs = list(
                .trainArgs=list(
                  formula=formula(Type~.),
                  cost=100)),
              metadata = list(
                reweighterInputPreds="preds",
                reweighterInputResponse="truth",
                reweighterInputWts="wts",
                reweighterOutputWts="wts")
              )


## ----, cache=FALSE-------------------------------------------------------
boostr::weightedAggregator


## ----, aggregatorExample, cache=FALSE------------------------------------
exoticAggr <- function(ensemble, estimator) {
  f <- ensemble[[estimator]]
  function(newdata) f(newdata)
}

boostr::boost(x = list(train = e1071::svm), B = 3,
              aggregator = exoticAggr,
              reweighter = boostr::arcfsReweighter,
              data = Glass,
              .procArgs = list(
                .trainArgs=list(
                  formula=formula(Type~.),
                  cost=100)),
              metadata = list(.inputEnsemble = "ensemble"),
              .boostBackendArgs = list(
                .aggregatorArgs = list(estimator = 2))
              )


## ----, cache=FALSE-------------------------------------------------------
boostr::defaultOOBPerformanceAnalysis


