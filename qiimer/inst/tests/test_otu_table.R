context('read_qiime_otu_table')

o <- read_qiime_otu_table(
  system.file("testdata", "otu_table.txt", package="qiimer"))

test_that("OTU IDs are parsed correctly", {
  expect_equal(o$otu_ids, c("0", "1", "2", "3", "5"))
})

test_that("Sample IDs are parsed correctly", {
  expect_equal(o$sample_ids, c("A.1", "A.2", "B.1", "B.2", "C.1"))
})

test_that("Counts are parsed correctly", {
  expect_equal(unname(o$counts[1,]), c(5, 10, 2, 0, 3))
})

test_that("Metadata is parsed correctly", {
  expect_equal(unname(o$metadata[1]), "Bacteria")
})

context('split_assignments')

a <- o$metadata

test_that('assignments are split into named ranks', {
  expect_equal(split_assignments(a)$Kingdom, factor(rep("Bacteria", 5)))
})

test_that('expected ranks are present', {
  expect_equal(colnames(split_assignments(a)), taxonomic_ranks[1:6])
})

test_that('assignment for OTU 3 is split correctly', {
  expect_equal(as.character(split_assignments(a)["3","Genus"]), "Weissella")
})

context('simplify_assignments')

sa <- split_assignments(a)

test_that('assignments before rank1 contain the full set of taxa', {
  expect_equal(
    simplify_assignments(sa, rank1="Class")["2"],
    structure("Bacteria Firmicutes", names="2"))
})

test_that('assignments up to rank1 contain only the taxon of rank1', {
  expect_equal(
    simplify_assignments(sa)["2"],
    structure("Firmicutes", names="2"))
})

test_that('assignments before rank2 contain rank1 plus most specific taxon', {
  expect_equal(
    simplify_assignments(sa)["5"],
    structure("Firmicutes Peptococcaceae", names="5"))
})

test_that('assignments past rank2 contain rank1 + rank2 taxa', {
  expect_equal(
    simplify_assignments(sa, rank2="Order")["5"],
    structure("Firmicutes Clostridiales", names="5"))
})

context('otu_heatmap')

ia <- simplify_assignments(sa)

test_that('otu_heatmap retains all taxa when threshold = 0', {
  expect_equal(nrow(otu_heatmap(o$counts, ia, plot=F)), 4)
})

test_that('otu_heatmap discards taxa below threshold', {
  expect_equal(nrow(otu_heatmap(o$counts, ia, threshold=5, plot=F)), 2)
})

test_that('otu_heatmap computes proportions before removing taxa', {
  h5 <- otu_heatmap(o$counts, ia, threshold=5, plot=F)
  h2 <- otu_heatmap(o$counts, ia, plot=F)
  expect_equal(h5, h2[rownames(h5),])
})

context('saturated_rainbow')

test_that('first element of saturated rainbow is white', {
  expect_equal(saturated_rainbow(20)[1], "#FFFFFFFF")
})

test_that('last element of saturated rainbow is red', {
  expect_equal(saturated_rainbow(20)[20], "#FF0000FF")
})

test_that('number of colors is equal to n', {
  expect_equal(length(saturated_rainbow(20)), 20)
})

