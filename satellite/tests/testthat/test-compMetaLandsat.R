context("Landsat metadata")

test_that("compFilePathLandsat works as expected for Landsat 7", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LE7*.tif"), 
                      full.names = TRUE)
  meta <- compFilePathLandsat(files)
  
  expect_equal(as.character(meta$BID[2]), "2")
  expect_equal(as.character(meta$BID[6]), "6_VCID_1")
  expect_equal(as.character(meta$BID[9]), "8")
  
  expect_equal(as.character(meta$BCDE[2]), "B002n")
  expect_equal(as.character(meta$BCDE[6]), "B0061")
  expect_equal(as.character(meta$BCDE[9]), "B008n")
  
  expect_equal(basename(as.character(meta$FILE[2])), 
               "LE71950252001211EDC00_B2.tif")
  expect_equal(basename(as.character(meta$FILE[6])), 
               "LE71950252001211EDC00_B6_VCID_1.tif")
  expect_equal(basename(as.character(meta$FILE[9])), 
               "LE71950252001211EDC00_B8.tif")
})

test_that("compFilePathLandsat works as expected for Landsat 8", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LC8*.tif"), 
                      full.names = TRUE)
  meta <- compFilePathLandsat(files)
  
  expect_equal(as.character(meta$BID[2]), "10")
  expect_equal(as.character(meta$BID[8]), "6")
  expect_equal(as.character(meta$BID[12]), "QA")
  
  expect_equal(as.character(meta$BCDE[2]), "B010n")
  expect_equal(as.character(meta$BCDE[8]), "B006n")
  expect_equal(as.character(meta$BCDE[12]), "B0QAn")
  
  expect_equal(basename(as.character(meta$FILE[2])), 
               "LC81950252013188LGN00_B10.tif")
  expect_equal(basename(as.character(meta$FILE[8])), 
               "LC81950252013188LGN00_B6.tif")
  expect_equal(basename(as.character(meta$FILE[12])), 
               "LC81950252013188LGN00_BQA.tif")
})


test_that("compMetaLandsat works as expected for Landsat 7", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LE7*.tif"), 
                      full.names = TRUE)
  meta <- compMetaLandsat(files)  
  
  expect_equal(as.character(meta$BID[2]), "2")
  expect_equal(as.character(meta$BID[6]), "6_VCID_1")
  expect_equal(as.character(meta$BID[9]), "8")
  
  expect_equal(as.character(meta$BCDE[2]), "B002n")
  expect_equal(as.character(meta$BCDE[6]), "B0061")
  expect_equal(as.character(meta$BCDE[9]), "B008n")
  
  expect_equal(basename(as.character(meta$FILE[2])), 
               "LE71950252001211EDC00_B2.tif")
  expect_equal(basename(as.character(meta$FILE[6])), 
               "LE71950252001211EDC00_B6_VCID_1.tif")
  expect_equal(basename(as.character(meta$FILE[9])), 
               "LE71950252001211EDC00_B8.tif")
})

test_that("compMetaLandsat works as expected for Landsat 8", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LC8*.tif"), 
                      full.names = TRUE)
  meta <- compMetaLandsat(files)  
  
  expect_equal(as.character(meta$BID[2]), "2")
  expect_equal(as.character(meta$BID[10]), "10")
  expect_equal(as.character(meta$BID[12]), "QA")
  
  expect_equal(as.character(meta$BCDE[2]), "B002n")
  expect_equal(as.character(meta$BCDE[10]), "B010n")
  expect_equal(as.character(meta$BCDE[12]), "B0QAn")
  
  expect_equal(basename(as.character(meta$FILE[2])), 
               "LC81950252013188LGN00_B2.tif")
  expect_equal(basename(as.character(meta$FILE[10])), 
               "LC81950252013188LGN00_B10.tif")
  expect_equal(basename(as.character(meta$FILE[12])), 
               "LC81950252013188LGN00_BQA.tif")
})

test_that("compMetaLandsat works as expected for Landsat 8 metafile", {
  path <- system.file("extdata",
                      package = "satellite")
  files <- list.files(path,
                      pattern = glob2rx("LC8*_MTL.txt"),
                      full.names = TRUE)
  meta <- compMetaLandsat(files)
  
  expect_equal(as.character(meta$BID[2]), "2")
  expect_equal(as.character(meta$BID[10]), "10")
  expect_equal(as.character(meta$BID[12]), "QA")
  
  expect_equal(as.character(meta$BCDE[2]), "B002n")
  expect_equal(as.character(meta$BCDE[10]), "B010n")
  expect_equal(as.character(meta$BCDE[12]), "B0QAn")
})
