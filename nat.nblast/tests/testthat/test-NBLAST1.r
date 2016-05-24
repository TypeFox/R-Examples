context("NBLAST v1")

testneurons <- readRDS('testdata/testneurons.rds')

test_that("nblast v1 produces expected scores", {
  scores <- nblast(testneurons[[1]], testneurons, version=1)
  scores.expected <- structure(c(0, 0.999980027716081, 0.999848861700461, 0.979269468123817, 0.999876876458676), .Names = c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001", "5HT1bMARCM-F000003_seg001", "5HT1bMARCM-F000004_seg001", "5HT1bMARCM-F000005_seg001"))
  expect_equal(scores, scores.expected)
})

test_that("nblast v1 with alpha produces expected scores", {
  scores <- nblast(testneurons[[1]], testneurons, version=1, UseAlpha=TRUE)
  scores.expected <- structure(c(0.202066124006236, 0.999982776382468, 0.999870496367803, 0.983132424575894, 0.999902562197814), .Names = c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001", "5HT1bMARCM-F000003_seg001", "5HT1bMARCM-F000004_seg001", "5HT1bMARCM-F000005_seg001"))
  expect_equal(scores, scores.expected)
})

test_that("nblast v1 handles a neuronlist as query", {
  scores <- nblast(testneurons[1:2], testneurons, version=1)
  scores.expected <- structure(c(0, 0.999980027716081, 0.999848861700461, 0.979269468123817, 0.999876876458676, 0.999947995458016, 0, 0.999999998533649, 0.999980102828588, 0.99999999937142), .Dim = c(5L, 2L), .Dimnames = list(c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001", "5HT1bMARCM-F000003_seg001", "5HT1bMARCM-F000004_seg001", "5HT1bMARCM-F000005_seg001"), c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001")))
  expect_equal(scores, scores.expected)
})

test_that("we can calculate normalised nblast v1 scores", {
  scores <- nblast(testneurons[1:2], testneurons, version=1)
  scores.norm=scale(scores, scale=c(scores[1,1],scores[2,2]), center = FALSE)
  expect_equivalent(nblast(testneurons[1:2], testneurons, version=1, normalised=TRUE), scores.norm)
})
