context("boost.list")

newPred <- function(obj, new) {
  predict(obj, new)
}

predMetadata <- c(modelName="obj", predictionSet="new")

B <- 3

test_that("boost.list works with homemade reweighters and aggregators", {
  
  metadata <- c(predMetadata, testReweighterMetadata, testAggregatorMetadata)

  boostedSVM <- boost(list(train=svm, predict=newPred),
                      B=B,
                      reweighter=testReweighter,
                      aggregator=testAggregator,
                      data=df,
                      .procArgs=testSVMProcArgs,
                      .boostBackendArgs = list(
                        .reweighterArgs=list(fakeStuff=77)),
                      metadata=metadata)
  
  expect_that(inherits(boostedSVM, 'boostr'), is_true())
  expect_that(length(ensembleEstimators(boostedSVM)), equals(B))
})

test_that("boost.list works with premade reweighters and aggregators", {
  metadata <- list( # prediction metadata
                   modelName="obj", predictionSet="new")

  ba <- list(
          arcfs = list(
            reweighter=arcfsReweighter,
            aggregator=arcfsAggregator,
            metadata=metadata),
          arcx4 = list(
            reweighter=arcx4Reweighter,
            aggregator=arcx4Aggregator,
            metadata=metadata,
            .boostBackendArgs=list(
              .reweighterArgs=list(m=0))),
          adaboost = list(
            reweighter=adaboostReweighter,
            metadata=metadata,
            aggregator=adaboostAggregator))
  
  lapply(ba, function(args) {
    out <- do.call(boost, c(list(x=list(train=svm, predict=newPred)), args,
                            B=B, list(.procArgs=testSVMProcArgs, data=df)))
    
    expect_that(inherits(out, 'boostr'), is_true())
    expect_that(length(ensembleEstimators(out)), equals(B))
  })
})

test_that("boost.list works with alternative performance analysis", {
  
  metadata <- c(predMetadata)
  
  boostedSVM <- boost(list(train=svm, predict=newPred),
                      B=B,
                      reweighter=vanillaBagger,
                      aggregator=vanillaAggregator,
                      data=df,
                      .procArgs=testSVMProcArgs,
                      analyzePerformance=testPerfAnalyzer,
                      .boostBackendArgs=list(
                        .analyzePerformanceArgs=list(Xx="77"),
                        .reweighterArgs=list(fakeStuff=77)),
                      metadata=metadata)
  
  expect_that(inherits(boostedSVM, 'boostr'), is_true())
  reweighterOutput <- attr(boostedSVM, "reweighterOutput")
  expect_that(length(ensembleEstimators(boostedSVM)), equals(B))
  expect_that(estimatorPerformance(boostedSVM),
              is_equivalent_to(matrix("77", nrow=B, ncol=1)))
  
  
  ### do this again with  testPerfAnalyzer2
  metadata <- c(predMetadata, testPerfAnalyzer2Metadata)
  
  boostedSVM <- boost(list(train=svm, predict=newPred),
                      B=B,
                      reweighter=vanillaBagger,
                      aggregator=vanillaAggregator,
                      data=df,
                      .procArgs=testSVMProcArgs,
                      analyzePerformance=testPerfAnalyzer2,
                      .boostBackendArgs=list(
                        .analyzePerformanceArgs=list(zeta="77"),
                        .reweighterArgs=list(fakeStuff=77)),
                      metadata=metadata)
  
  expect_that(inherits(boostedSVM, 'boostr'), is_true())
  reweighterOutput <- attr(boostedSVM, "reweighterOutput")
  expect_that(length(ensembleEstimators(boostedSVM)), equals(B))
  expect_that(length(estimatorPerformance(boostedSVM)), equals(B))
  expect_that(estimatorPerformance(boostedSVM)[[1]]$z, equals("77"))
})