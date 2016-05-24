context("boost.function")

B <- 3

test_that("vignette examples work", {
  library(mlbench)
  data(Glass)
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
  
  out <- boostWithArcX4(x = kNN_EstimationProcedure,
                         B = 3,
                         data = Glass,
                         metadata = list(learningSet="learningSet"),
                         .procArgs = list(k=5),
                         .boostBackendArgs = list(
                           .subsetFormula=formula(Type~.))
                        ) 
})

test_that("boost.function works with homemade reweighters and aggregators", {
  metadata <- c(testkNNProcMetadata, testReweighterMetadata,
                testAggregatorMetadata)
  
  boostedkNN <- boost(testkNNProc,
                      B=B,
                      reweighter=testReweighter,
                      aggregator=testAggregator,
                      data=df,
                      .procArgs=testKNNProcArgs,
                      .boostBackendArgs = list(
                        .reweighterArgs=list(fakeStuff=77)),
                      metadata=metadata)
  
  expect_that(inherits(boostedkNN, 'boostr'), is_true())
  expect_that(length(ensembleEstimators(boostedkNN)), equals(B))
})

test_that("boost.function works with premade reweighters and aggregators", {
  metadata <- c(testkNNProcMetadata)

  ba <- list(arcfs = list(
              reweighter=arcfsReweighter,
              aggregator=arcfsAggregator,
              metadata=metadata),
             arcx4 = list(
               reweighter = arcx4Reweighter,
               aggregator = arcx4Aggregator,
               metadata = metadata,
               .boostBackendArgs = list(.reweighterArgs=list(m=0))),
             adaboost = list(
               reweighter=adaboostReweighter,
               metadata=metadata,
               aggregator=adaboostAggregator))
  
  lapply(ba, function(args) {

    out <- do.call(boost, c(x=testkNNProc, args, B=B, 
                            list(.procArgs=testKNNProcArgs, data=df)))
    
    expect_that(inherits(out, 'boostr'), is_true())
    expect_that(length(ensembleEstimators(out)), equals(B))
  })
})

test_that("boost.function works alternative performance analysis", {
  metadata <- c(testkNNProcMetadata)
  
  boostedkNN <- boost(testkNNProc,
                      B=B,
                      reweighter=vanillaBagger,
                      aggregator=vanillaAggregator,
                      data=df,
                      .procArgs=testKNNProcArgs,
                      metadata=metadata,
                      .boostBackendArgs=list(
                        .analyzePerformanceArgs = list(Xx="77"),
                        .reweighterArgs=list(fakeStuff=77)),
                      analyzePerformance=testPerfAnalyzer)
  
  expect_that(inherits(boostedkNN, 'boostr'), is_true())
  reweighterOutput <- attr(boostedkNN,  "reweighterOutput")
  expect_that(length(ensembleEstimators(boostedkNN)), equals(B))
  expect_that(estimatorPerformance(boostedkNN),
              is_equivalent_to(matrix("77", nrow=B, ncol=1)))
  
  ### do it again with testPerfAnalyzer2
  
  metadata <- c(metadata, testPerfAnalyzer2Metadata)
  
  boostedkNN <- boost(testkNNProc,
                      B=B,
                      data=df,
                      reweighter=vanillaBagger,
                      aggregator=vanillaAggregator,
                      .procArgs=testKNNProcArgs,
                      analyzePerformance=testPerfAnalyzer2,
                      .boostBackendArgs=list(
                        .analyzePerformanceArgs = list(zeta="77"),
                        .reweighterArgs = list(fakeStuff=77)),
                      metadata=metadata)
  
  expect_that(inherits(boostedkNN, 'boostr'), is_true())
  reweighterOutput <- attr(boostedkNN,  "reweighterOutput")
  expect_that(length(ensembleEstimators(boostedkNN)), equals(B))
  expect_that(length(estimatorPerformance(boostedkNN)), equals(B))
  expect_that(estimatorPerformance(boostedkNN)[[1]]$z, equals("77"))
})