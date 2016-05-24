# devtools::test(".", "satellite")
context("satellite")

test_that("satellite works as expected for Landsat 7  files", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LE7*.tif"), 
                      full.names = TRUE)
  
  sat <- satellite(files)

  expect_equal(basename(attr(getSatDataLayers(sat)[[1]], "file")@name), "LE71950252001211EDC00_B1.tif")
  expect_equal(basename(attr(getSatDataLayers(sat)[[2]], "file")@name), "LE71950252001211EDC00_B2.tif")
  expect_equal(basename(attr(getSatDataLayers(sat)[[6]], "file")@name), "LE71950252001211EDC00_B6_VCID_1.tif")
  expect_equal(basename(attr(getSatDataLayers(sat)[[9]], "file")@name), "LE71950252001211EDC00_B8.tif")
  
  expect_equal(as.character(getSatMeta(sat)$BID[[1]]), "1")
  expect_equal(as.character(getSatMeta(sat)$BID[[2]]), "2")
  expect_equal(as.character(getSatMeta(sat)$BID[[6]]), "6_VCID_1")
  expect_equal(as.character(getSatMeta(sat)$BID[[9]]), "8")
  
  expect_equal(as.character(getSatMeta(sat)$BCDE[[1]]), "B001n")
  expect_equal(as.character(getSatMeta(sat)$BCDE[[2]]), "B002n")
  expect_equal(as.character(getSatMeta(sat)$BCDE[[6]]), "B0061")
  expect_equal(as.character(getSatMeta(sat)$BCDE[[9]]), "B008n")
})

test_that("satellite works as expected for Landsat 8 files", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LC8*.tif"), 
                      full.names = TRUE)
  sat <- satellite(files)
  
  expect_equal(names(getSatDataLayers(sat)[[1]]), "B001n")
  expect_equal(names(getSatDataLayers(sat)[[2]]), "B002n")
  expect_equal(names(getSatDataLayers(sat)[[8]]), "B008n")
  expect_equal(names(getSatDataLayers(sat)[[10]]), "B010n")
  expect_equal(names(getSatDataLayers(sat)[[12]]), "B0QAn")
  
  expect_equal(as.character(getSatMeta(sat)$BID[[1]]), "1")
  expect_equal(as.character(getSatMeta(sat)$BID[[2]]), "2")
  expect_equal(as.character(getSatMeta(sat)$BID[[6]]), "6")
  expect_equal(as.character(getSatMeta(sat)$BID[[12]]), "QA")
  
  expect_equal(as.character(getSatMeta(sat)$BCDE[[1]]), "B001n")
  expect_equal(as.character(getSatMeta(sat)$BCDE[[2]]), "B002n")
  expect_equal(as.character(getSatMeta(sat)$BCDE[[6]]), "B006n")
  expect_equal(as.character(getSatMeta(sat)$BCDE[[12]]), "B0QAn")
})


test_that("satellite works as expected for stacks", {
  sat <- satellite(l8)
})