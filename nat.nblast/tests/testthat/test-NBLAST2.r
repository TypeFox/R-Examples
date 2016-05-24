context("NBLAST v2")

testneurons <- readRDS('testdata/testneurons.rds')

test_that("nblast v2 produces expected scores", {
  scores <- nblast(testneurons[[1]], testneurons, version=2)
  scores.expected <- structure(c(43518.2468824115, -34509.2778341796, -36608.1290149345, -18366.0246697251, -36782.7167494325), .Names = c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001", "5HT1bMARCM-F000003_seg001", "5HT1bMARCM-F000004_seg001", "5HT1bMARCM-F000005_seg001"))

  expect_equal(scores, scores.expected)
})

test_that("nblast v2 with alpha produces expected scores", {
  scores <- nblast(testneurons[[1]], testneurons, version=2, UseAlpha=TRUE)
  scores.expected <- structure(c(22393.1456610767, -33103.620929615, -35603.5457297204, -18717.4612913414, -35799.9999942706), .Names = c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001", "5HT1bMARCM-F000003_seg001", "5HT1bMARCM-F000004_seg001", "5HT1bMARCM-F000005_seg001"))
  expect_equal(scores, scores.expected)
})

test_that("nblast v2 handles a neuronlist as query", {
  scores <- nblast(testneurons[1:2], testneurons, version=2)
  scores.expected <- structure(c(43518.2468824115, -34509.2778341796, -36608.1290149345, -18366.0246697251, -36782.7167494325, -14764.1401323446, 30648.4172626457, -22373.118229109, -22772.3018332237, -22918.5324546932), .Dim = c(5L, 2L), .Dimnames = list(c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001", "5HT1bMARCM-F000003_seg001", "5HT1bMARCM-F000004_seg001", "5HT1bMARCM-F000005_seg001"), c("5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001")))
  expect_equal(scores, scores.expected)
})

test_that("we can calculate normalised nblast v2 scores", {
  scores <- nblast(testneurons[1:2], testneurons, version=2)
  scores.norm=scale(scores, scale=c(scores[1,1],scores[2,2]), center = FALSE)
  expect_equivalent(nblast(testneurons[1:2], testneurons, normalised=TRUE), scores.norm)
})

test_that("we can calculate scores for regular neurons",{
  library(nat)
  nblast_allbyall(Cell07PNs[1:4])
})

test_that("we can calculate scores using getOption('nat.default.neuronlist')", {
  op=options(nat.default.neuronlist='testneurons')
  on.exit(options(op))
  scores <- nblast(testneurons[1:2], testneurons, version=2)
  expect_equal(nblast(testneurons[1:2], version=2), scores)

  options(nat.default.neuronlist='testneuronsrhubarb')
  expect_error(nblast(testneurons[1:2], version=2))
})

test_that("we can handle OmitFailures", {
  testneurons.err=testneurons
  testneurons.err[[2]]='rhubarb'
  expect_error(nblast(testneurons.err[1:2], testneurons))
  scores <- nblast(testneurons[1:2], testneurons, version=2)
  scores.err=scores
  scores.err[,2]=NA_real_
  expect_equal(nblast(testneurons.err[1:2], testneurons, OmitFailures = FALSE),
               scores.err)
  expect_error(nblast(testneurons.err[1:2], testneurons, OmitFailures = TRUE),
               "not yet implemented")

})

test_that("we can handle all combinations of dotprops and neurons, both as neuronlists and singly", {
  nblast(testneurons[[1]], Cell07PNs[1:3])
  nblast(testneurons[1:3], Cell07PNs[1:3])
  nblast(testneurons[[1]], Cell07PNs[[1]])
  nblast(testneurons[1:3], Cell07PNs[[1]])

  nblast(Cell07PNs[[1]], testneurons[1:3])
  nblast(Cell07PNs[1:3], testneurons[1:3])
  nblast(Cell07PNs[[1]], testneurons[[1]])
  nblast(Cell07PNs[1:3], testneurons[[1]])
})
