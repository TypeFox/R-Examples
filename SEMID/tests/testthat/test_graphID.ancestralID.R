library(SEMID)
context("Testing known examples on graphID.ancestralID.")

source("graphExamples.R")

test_that("graphID.ancestralID returns correct value for known examples.", {
  for (i in 1:length(graphExamples)) {
    graphExample = graphExamples[[i]]
    L = graphExample$L
    O = graphExample$O
    htcId = graphExample$htcId
    ancId = graphExample$ancId
    m = nrow(L)
    if (!is.dag(graph.adjacency(L))) {
      next
    }

    result = graphID.ancestralID(L, O)

    if (!is.null(ancId)) {
      if (ancId == 1) {
        expect_equal(sort(result), 1:m)
      } else {
        expect_true(all(result %in% 1:m) && length(result) < m)
      }
    }

    if (htcId == 0) {
      expect_true(all(result %in% 1:m) && length(result) < m)
    } else if (htcId == 1) {
      expect_equal(sort(result), 1:m)
    } else {
      expect_true(all(graphID.htcID(L, O) %in% result))
    }
  }
})
