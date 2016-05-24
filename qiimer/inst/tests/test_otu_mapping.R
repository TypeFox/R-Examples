context('read_qiime_otu_mapping')

m <- read_qiime_otu_mapping(
  system.file("testdata", "otu_mapping.txt", package="qiimer"))

test_that("OTU IDs are parsed correctly", {
  expect_equal(names(m), c("0", "1", "2", "3", "5"))
})

test_that("Sequence IDs are parsed correctly", {
  expect_equal(m[["0"]], "A.1_216513")
})

test_that("Number of sequences is correct for each OTU", {
  expect_equal(
    sapply(m, length), 
    structure(c(1, 1, 1, 6, 5), names=c("0", "1", "2", "3", "5")))
})

context('make_otu_table')

expected_table <- as.table(matrix(
  data=c(
    1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 
    0, 0, 1, 0, 0, 
    6, 0, 0, 0, 0, 
    1, 1, 1, 1, 1), 
  nrow=5, 
  byrow=T,
  dimnames=list(
    OTU=c("0", "1", "2", "3", "5"),
    SampleID=c("A.1", "A.2", "B.1", "B.2", "C.1"))))

expect_equal(make_otu_table(m), expected_table)
