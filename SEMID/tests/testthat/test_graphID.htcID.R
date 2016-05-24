library(SEMID)
context("Testing that half-trek criterion function for generic identifiability works properly.")

source("graphExamples.R")

test_that("graphID.htcID returns correct value for known examples.", {
  for (i in 1:length(graphExamples)) {
    graphExample = graphExamples[[i]]
    L = graphExample$L
    O = graphExample$O
    htcId = graphExample$htcId
    m = nrow(L)

    result = graphID.htcID(L, O)
    if (htcId == 1) {
      expect_equal(sort(result), 1:m)
    } else {
      expect_true(length(result) == 0 ||
                    (all(result %in% 1:m) && length(result) != m))
    }
  }
})