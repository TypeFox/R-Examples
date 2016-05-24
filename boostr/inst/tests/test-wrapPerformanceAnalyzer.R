context("wrapPerformanceAnalyzer")

wrapTestPA2 <- do.call(wrapPerformanceAnalyzer,
                       c(analyzePerformance=testPerfAnalyzer2,
                         testPerfAnalyzer2Metadata))

test_that("wrapped performance analyzer has the right class", {
  expect_that(wrapTestPA2, is_a("performanceAnalyzer"))
})

metadata <- testPerfAnalyzer2Metadata
test_that("wrapPerformanceAnalyzer is compatible with boostBackend", {
  phi <- boostBackend(B=3,
                      initialWeights=rep.int(1, nrow(df)),
                      reweighter=vanillaBagger,
                      aggregator=vanillaAggregator,
                      proc=testSVM,
                      .procArgs=testSVMProcArgs,
                      data=df,
                      analyzePerformance=wrapTestPA2,
                      .analyzePerformanceArgs = list(zeta=0))
  
  
  lapply(estimatorPerformance(phi), function(x) {
    expect_that(x$z, equals(0))
  })                   
})

test_that("wrapPerformanceAnalyzer is compatible inside boost", {
  phi <- boost(x=testSVM, B=3,
               reweighter=vanillaBagger,
               aggregator=vanillaAggregator,
               data=df,
               .procArgs=testSVMProcArgs,
               analyzePerformance=wrapTestPA2,
               .boostBackendArgs = list(
                 .analyzePerformanceArgs=list(zeta=0)),
               metadata=metadata)

  
  lapply(estimatorPerformance(phi), function(x) {
    expect_that(x$z, equals(0))
  })  
})

test_that("wrapPerformanceAnalyzer gets called from inside boost", {
  phi <- boost(x=testSVM,
               B=3,
               reweighter=vanillaBagger,
               aggregator=vanillaAggregator,
               data=df,
               .procArgs=testSVMProcArgs,
               analyzePerformance=testPerfAnalyzer2,
               .boostBackendArgs = list(
                 .analyzePerformanceArgs=list(zeta=0)),
               metadata=metadata)
  
  lapply(estimatorPerformance(phi), function(x) {
    expect_that(x$z, equals(0))
  })  
})