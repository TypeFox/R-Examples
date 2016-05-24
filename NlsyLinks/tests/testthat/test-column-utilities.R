
# LoadPairFile <- function( ) {
#   #directory <- "F:/Projects/Nls/Links2011/Analysis/Df/2012-01-13"
#   #pathLinks <- file.path(directory, "Links2011V28.csv")
#   #dsLinks <- read.csv(pathLinks)
#   return( Links79Pair )
# }
fileNameGen2 <- "gen2-birth.csv"

# print(basename(normalizePath(".")))
{
  if( basename(normalizePath("."))=="NlsyLinks" ) {
    directoryForExpectedVectors <- "./inst/tests"
  }
  else if( basename(normalizePath("."))=="testthat" ) {
    directoryForExpectedVectors <- "."
  }
  else {
    stop("The working directory is not recognized by this test fixture.")
  }
}
source(file.path(directoryForExpectedVectors, "expected-vectors.R"))


###########
context("CreateSubjectTag")
###########
test_that("CreateSubjectTag -Scenario 1", {
  ids <- c(1:10, 1:10)
  generation <- c(rep(1, 10), rep(2, 10))
  expected <- c(1:10*100, 1:10)
  expect_equal(expected, CreateSubjectTag(ids, generation))
  expect_equivalent(expected, CreateSubjectTag(ids, generation))
})
test_that("CreateSubjectTag -Scenario 2", {
  ids <- c(71:80, 1:10)
  generation <- c(rep(2, 10), rep(1, 10))
  expected <- c(71:80, 1:10*100)
  expect_equal(expected, CreateSubjectTag(ids, generation))
  expect_equivalent(expected, CreateSubjectTag(ids, generation))
})
test_that("CreateSubjectTag -Scenario 3", {
  ids <- c(NA, NA, 71:80, NA, 1:10, NA)
  generation <- c(rep(2, 12), rep(1, 12))
  expected <- c(NA, NA,71:80, NA, 1:10*100, NA)
  expect_equal(expected, CreateSubjectTag(ids, generation))
  expect_equivalent(expected, CreateSubjectTag(ids, generation))
})
test_that("CreateSubjectTag -Scenario 4", {
  ids <- c(71:82, 1:12)
  generation <- c(NA, NA, rep(2, 10), NA, rep(1, 10), NA)
  expected <- c(NA, NA,73:82, NA, 2:11*100, NA)
  expect_equal(expected, CreateSubjectTag(ids, generation))
  expect_equivalent(expected, CreateSubjectTag(ids, generation))
})
test_that("CreateSubjectTag -Scenario 5", {
  ids <- c(71:82, 10001:10012)
  generation <- c(NA, NA, rep(1, 10), NA, rep(2, 10), NA)
  expected <- c(NA, NA,73:82*100, NA, 10002:10011, NA)
  expect_equal(expected, CreateSubjectTag(ids, generation))
  expect_equivalent(expected, CreateSubjectTag(ids, generation))
})

test_that("CreateSubjectTag -Scenario 6", {
  ids <- c(71:82, 10001:10012)
  generation <- c( rep(1, 12),rep(2, 12))
  expected <- c(71:82*100, 10001:10012)
  expect_equal(expected, CreateSubjectTag(ids, generation))
  expect_equivalent(expected, CreateSubjectTag(ids, generation))
})

test_that("CreateSubjectTag -With ExtraOutcomes79", {  
  actual <- CreateSubjectTag(subjectID=ExtraOutcomes79$SubjectID, generation=ExtraOutcomes79$Generation)
  expected <- ExpectedSubjectTags
  
  expect_equal(expected, actual)
  #The PrintVector function is in "./Content/DeveloperUtilities.R"
  #cat(PrintVector(ExtraOutcomes79$SubjectTag))  
})




###########
context("ExtractColumnExists")
###########
test_that("Nlsy79Gen2", {
  #filePathGen2 <- file.path(path.package("NlsyLinks"), "extdata", "gen2-birth.csv") #"./Datasets/gen2-birth.csv"  
  filePathGen2 <- file.path(devtools::inst("NlsyLinks"), "extdata", "gen2-birth.csv") #"./Datasets/gen2-birth.csv"  
  expectedColumNames <- c("C0000100", "C0000200", "C0005300", "C0005400", "C0005700", "C0328000", "C0328600", "C0328800")
  ds <- read.csv(filePathGen2)
  expectedIndex <- 0
  for( expectedColumnName in expectedColumNames ) {
    expectedIndex <- expectedIndex + 1
    expect_equal(VerifyColumnExists(ds, expectedColumnName), expected=expectedIndex, info=paste("The column '", expectedColumnName, "' should be found."))
  }
})

###########
context("Rename Nlsy Column")
###########
test_that("RenameNlsyColumn", {
  filePathGen2 <- file.path(devtools::inst("NlsyLinks"), "extdata", "gen2-birth.csv") #"./Datasets/gen2-birth.csv"  
  ds <- read.csv(filePathGen2)
  originalColumNames <- c("C0000100", "C0000200", "C0005300", "C0005400", "C0005700", "C0328000", "C0328600", "C0328800")
  newColumnNames <- c("SubjectID", "MotherID", "Race", "Gender", "Yob", "GestationWeeks", "BirthWeightInOunces", "BirthLengthInInches")
  for( columnIndex in seq_along(originalColumNames) ) {
    ds <- RenameNlsyColumn(dataFrame=ds, nlsyRNumber=originalColumNames[columnIndex], newColumnName=newColumnNames[columnIndex])
    #expect_equal(VerifyColumnExists(ds, expectedColumnName), expected=expectedIndex, paste("The column '", expectedColumnName, "' should be found."))
  }
  expect_equal(colnames(ds), expected=newColumnNames, info="The renamed columns should be correct.")
  
})


# 
# test_that("Zero rows", {
#   dsLinks <- LoadPairFile()
#   dsLinks <- dsLinks[0,]
#   expect_error(ValidatePairLinks(dsLinks), "The linksPair file should have at least one row, but does not.")
# })
# 
# test_that("Bad SubjectTag_S1", {
#   dsLinks <- LoadPairFile()
#   expect_true(ValidatePairLinks(dsLinks))
#   colnames(dsLinks)[colnames(dsLinks)=="SubjectTag_S1"] <- "Bad"
#   expect_error(ValidatePairLinks(dsLinks), "The column 'SubjectTag_S1' should exist in the linksPair file, but does not.")
# })
# 
# test_that("Bad SubjectTag_S2", {
#   dsLinks <- LoadPairFile()
#   expect_true(ValidatePairLinks(dsLinks))
#   colnames(dsLinks)[colnames(dsLinks)=="SubjectTag_S2"] <- "Bad"
#   expect_error(ValidatePairLinks(dsLinks), "The column 'SubjectTag_S2' should exist in the linksPair file, but does not.")
# })
# 
# test_that("Bad R", {
#   dsLinks <- LoadPairFile()
#   expect_true(ValidatePairLinks(dsLinks))
#   colnames(dsLinks)[colnames(dsLinks)=="R"] <- "Bad"
#   expect_error(ValidatePairLinks(dsLinks), "The column 'R' should exist in the linksPair file, but does not.")
# })
# 
# # test_that("Bad MultipleBirth", {
# #   dsLinks <- LoadPairFile()
# #   expect_true(ValidatePairLinks(dsLinks))
# #   colnames(dsLinks)[colnames(dsLinks)=="MultipleBirth"] <- "Bad"
# #   expect_error(ValidatePairLinks(dsLinks), "The column 'MultipleBirth' should exist in the linksPair file, but does not.")
# # })
