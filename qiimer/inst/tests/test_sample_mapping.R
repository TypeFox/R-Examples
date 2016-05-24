context('read_qiime_mapping_file')

s <- read_qiime_mapping_file(
  system.file("testdata", "sample_map.txt", package="qiimer"))

test_that("Column names are parsed correctly", {
  expect_equal(colnames(s), c(
    "SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Diet", 
    "Description"))
})

test_that("Sample IDs are parsed correctly", {
  expect_equal(class(s$SampleID), "character")
  expect_equal(s$SampleID, c("A.1", "A.2", "B.1", "B.2", "C.1"))
})
