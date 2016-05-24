library(devtools)
# filePathGen2 <- "./../extdata/gen2-birth.csv" #"./Datasets/gen2-birth.csv"
# fileNameGen2 <- "gen2-birth.csv"

# print(basename(normalizePath(".")))
# {
# if( basename(normalizePath("."))=="NlsyLinks" ) 
#   directory <- "./inst/extdata"
# else if( basename(normalizePath("."))=="tests" ) 
#   directory <- "./../extdata/"
# else
#   stop("The working directory is not recognized by this test fixture.")
# }
# basename(dirname(normalizePath(".")))
# basename((normalizePath(".")))

###########
context("Read CSV")
###########
test_that("Nlsy79Gen1Path", {  
  filePathGen1 <- file.path(devtools::inst("NlsyLinks"), "extdata", "gen1-life-course.csv")
#   dsExtract <- read.csv(filePathGen1)
  ds <- ReadCsvNlsy79Gen1(filePath=filePathGen1)
  
  expect_equal(object=min(ds$SubjectTag), expected=100, scale=1)
  expect_equal(object=max(ds$SubjectTag), expected=1268600, scale=1)
  expect_true(all(ds$Generation==1))
  expect_equal(object=nrow(ds), expected=12686, scale=1)
  expect_equal(object=ncol(ds), expected=13, scale=1)
})
test_that("Nlsy79Gen1DataFrame", {  
  filePathGen1 <- file.path(devtools::inst("NlsyLinks"), "extdata", "gen1-life-course.csv")
  dsRaw <- read.csv(filePathGen1)
  ds <- ReadCsvNlsy79Gen1(dsExtract=dsRaw)
  
  expect_equal(object=min(ds$SubjectTag), expected=100, scale=1)
  expect_equal(object=max(ds$SubjectTag), expected=1268600, scale=1)
  expect_true(all(ds$Generation==1))
  expect_equal(object=nrow(ds), expected=12686, scale=1)
  expect_equal(object=ncol(ds), expected=13, scale=1)
})

test_that("Nlsy79Gen2Path", {  
  #   ds <- ReadCsvNlsy79Gen2(filePath=file.path(directory, fileNameGen2))
  filePathGen2 <- file.path(devtools::inst("NlsyLinks"), "extdata", "gen2-birth.csv")
  ds <- ReadCsvNlsy79Gen2(filePath=filePathGen2)
  
  expect_equal(object=min(ds$SubjectTag), expected=201, scale=1)
  expect_equal(object=max(ds$SubjectTag), expected=1267501, scale=1)
  expect_true(all(ds$Generation==2))
  expect_equal(object=min(ds$SubjectTagOfMother), expected=200, scale=1)
  expect_equal(object=max(ds$SubjectTagOfMother), expected=1267500, scale=1)
  expect_equal(object=nrow(ds), expected=11495, scale=1)
  expect_equal(object=ncol(ds), expected=11, scale=1)
})

test_that("Nlsy79Gen2DataFrame", {  
  filePathGen2 <- file.path(devtools::inst("NlsyLinks"), "extdata", "gen2-birth.csv")
  dsRaw <- read.csv(filePathGen2)
  ds <- ReadCsvNlsy79Gen2(dsExtract=dsRaw)
  
  expect_equal(object=min(ds$SubjectTag), expected=201, scale=1)
  expect_equal(object=max(ds$SubjectTag), expected=1267501, scale=1)
  expect_true(all(ds$Generation==2))
  expect_equal(object=min(ds$SubjectTagOfMother), expected=200, scale=1)
  expect_equal(object=max(ds$SubjectTagOfMother), expected=1267500, scale=1)
  expect_equal(object=nrow(ds), expected=11495, scale=1)
  expect_equal(object=ncol(ds), expected=11, scale=1)
})
