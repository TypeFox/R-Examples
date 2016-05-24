# Example metadata for BIOM object
b_otu_ids <- c("denovo0", "denovo1", "denovo2", "denovo3", "denovo5")
b_sample_ids <- c("A.1", "A.2", "B.1", "B.2", "C.1")
b_taxonomy <- list(
  "Bacteria", 
  "Bacteria", 
  c("Bacteria", "Firmicutes"), 
  c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", 
    "Leuconostocaceae", "Weissella"), 
  c("Bacteria", "Firmicutes", "Clostridia", "Clostridiales", 
    "Peptococcaceae"))

# Prototype BIOM object without data
b <- list(type = "OTU table", shape=c(5, 5))
b$rows <- mapply(
  function (x, y) list(id=x, metadata=list(taxonomy=y)),
  b_otu_ids, b_taxonomy, SIMPLIFY=FALSE, USE.NAMES=FALSE)
b$columns <- lapply(b_sample_ids, function (x) list(id=x, metadata=NULL))

# BIOM object in dense format
db <- b
db$matrix_type <- "dense"
db$data <- list(
  c(5, 10, 2, 0, 3), 
  c(76, 23, 28, 43, 56), 
  c(2, 1, 0, 0, 0), 
  c(637, 61, 63, 77, 34), 
  c(0, 0, 0, 1, 0))

# BIOM object in sparse format
sb <- b
sb$matrix_type <- "sparse"
sb$data <- list(
  c(0, 0, 5), c(0, 1, 10), c(0, 2, 2), c(0, 4, 3),
  c(1, 0, 76), c(1, 1, 23), c(1, 2, 28), c(1, 3, 43), c(1, 4, 56),
  c(2, 0, 2), c(2, 1, 1),
  c(3, 0, 637), c(3, 1, 61), c(3, 2, 63), c(3, 3, 77), c(3, 4, 34),
  c(4, 3, 1))


context('biom_raw_data, dense format')

test_that("OTU IDs are parsed correctly", {
  expect_equal(rownames(biom_raw_data(db)), b_otu_ids)
})

test_that("Sequence IDs are parsed correctly", {
  expect_equal(colnames(biom_raw_data(db)), b_sample_ids)
})

test_that("OTU counts are correct", {
  expect_equal(
    biom_raw_data(db)["denovo0",], 
    structure(c(5, 10, 2, 0, 3), names=b_sample_ids))
  expect_equal(
    biom_raw_data(db)[,"A.1"], 
    structure(c(5, 76, 2, 637, 0), names=b_otu_ids))
})


context('biom_raw_data, sparse format')

test_that("OTU IDs are present", {
  expect_equal(levels(biom_raw_data(sb)$OTU), b_otu_ids)
})

test_that("Sequence IDs are present", {
  expect_equal(levels(biom_raw_data(sb)$SampleID), b_sample_ids)
})

test_that("OTU counts are in the same order as the biom file", {
  expect_equal(
    biom_raw_data(sb)$value, 
    c(5, 10, 2, 3, 76, 23, 28, 43, 56, 2, 1, 637, 61, 63, 77, 34, 1))
})

test_that("Values are correct if OTUs/samples do not appear in data", {
  ssb <- sb
  ssb$data <- ssb$data[8:11]
  expect_equal(
    biom_raw_data(ssb), 
    data.frame(
      OTU=c("denovo1", "denovo1", "denovo2", "denovo2"),
      SampleID=c("B.2", "C.1", "A.1", "A.2"),
      value=c(43, 56, 2, 1)))
})


context('biom_taxonomy')

test_that("Taxonomy is correct", {
  expect_equal(biom_taxonomy(sb)[[1]], "Bacteria")
})

test_that("Taxonomy list is labeled with OTU IDs", {
  expect_equal(biom_taxonomy(sb)[["denovo2"]], c("Bacteria", "Firmicutes"))
})
