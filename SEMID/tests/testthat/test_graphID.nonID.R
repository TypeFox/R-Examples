library(SEMID)
context("Testing that half-trek criterion function for generic non-identifiability works properly.")

source("graphExamples.R")

test_that("graphID.nonHtcID returns correct value for known examples.", {
  for (i in 1:length(graphExamples)) {
    graphExample = graphExamples[[i]]
    L = graphExample$L
    O = graphExample$O
    htcId = graphExample$htcId
    m = nrow(L)

    result = graphID.nonHtcID(L, O)
    if (htcId == 0) {
      expect_true(result)
    } else {
      expect_false(result)
    }
  }
})
