context("wrapAggregator")

wrappedAggregator <- do.call(wrapAggregator,
                             list(aggregator=testAggregator,testAggregatorMetadata))

phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(df)),
              reweighter=arcx4Reweighter,
              aggregator=arcx4Aggregator,
              proc=testSVM,
              .procArgs=testSVMProcArgs,
              data=df,
              .reweighterArgs=list(m=0))

ensemble <- ensembleEstimators(phi)

test_that("aggregator wrapper doesn't affect results", {
  set.seed(1234)
  originalEstimator <- testAggregator(ensemble)
  originalPreds <- originalEstimator(iris[1:15, ])
  
  set.seed(1234)
  wrappedEstimator <- wrappedAggregator(ensemble)
  wrappedPreds <- wrappedEstimator(iris[1:15, ])
  
  expect_that(originalPreds, equals(wrappedPreds))
})

test_that("wrapAggregator is compatible with boostBackend", {
  phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(df)),
                    reweighter=arcx4Reweighter,
                    aggregator=wrappedAggregator,
                    proc=testSVM,
                    .procArgs=testSVMProcArgs,
                    data=df,
                    .reweighterArgs=list(m=0))
  
  expect_that(is.null(phi), is_false())
})