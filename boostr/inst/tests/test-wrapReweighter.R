context("wrapReweighter")

set.seed(1234)
out1 <- testReweighter(preds, response, wts)

wrappedTestReweighter <- do.call(wrapReweighter,
                              c(reweighter=testReweighter, testReweighterMetadata))

set.seed(1234)
out2 <- wrappedTestReweighter(prediction=preds, response=response, weights=wts)

test_that("reweighter wrapper doesn't affect results", {
  sapply(seq.int(length(out1)), function(entry) {
    expect_that(out1[[entry]], is_equivalent_to(out2[[entry]]))
  })
})

test_that("wrapReweighter is compatible with boostBackend", {
  phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(df)),
                       reweighter=wrappedTestReweighter,
                       aggregator=arcx4Aggregator,
                       proc=testSVM,
                       .procArgs=testSVMProcArgs,
                       .reweighterArgs=list(m=0),
                       data=df,
                       .formatData=TRUE)
  
  expect_that(is.null(phi), is_false())
})