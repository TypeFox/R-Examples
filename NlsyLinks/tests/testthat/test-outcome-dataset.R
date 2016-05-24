
# LoadPairFile <- function( ) {
#   return( Links79Pair )
# }
LoadOutcomeFile <- function( ) {
  return( ExtraOutcomes79 )
}
LoadDefaultOutcomeNames <- function( ) {
  return( c("AfqtRescaled2006Gaussified", "HeightZGenderAge", 
            "WeightZGenderAge", "Afi", "Afm", "MathStandardized") )
}



###########
context("Validate Outcome Dataset")
###########
test_that("Normal Scenario", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)
  expect_true(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()))
})

test_that("Zero rows", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)  
  ds <- ds[0,]
  expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()), "The dsOutcome data frame should have at least one row, but does not.")
})

test_that("dsOutcome -Missing ", {
  ds <- LoadOutcomeFile()
  expect_error(ValidateOutcomeDataset( outcomeNames=LoadDefaultOutcomeNames()), "The parameter for 'dsOutcome' should be passed, but was not.")
})

# test_that("SubjectTag -Missing ", {
#   ds <- LoadOutcomeFile()
#   expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()), "The column 'SubjectTag' should exist in the data frame, but does not.")
# })

test_that("SubjectTag -bad mode", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- sample(letters, nrow(ds), replace=T)
  expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()), "The column 'SubjectTag' should have a 'numeric' mode, but does not.")
})

test_that("SubjectTag -negative values", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)
  ds$SubjectTag[33] <-  -ds$SubjectTag[33]
  expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()), "The column 'SubjectTag' should contain only positive values, but does not.")
})

test_that("SubjectTag -NA values", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)
  ds$SubjectTag[333] <- NA
  expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()), "The column 'SubjectTag' should not contain any NA values, but it does.")
})

test_that("SubjectTag -duplicate values -scenario 1", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)
  ds$SubjectTag[2] <- ds$SubjectTag[1]
  expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()), "The column 'SubjectTag' should not contain duplicated, but it does.")
})

test_that("SubjectTag -duplicate values -scenario 2", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)
  ds$SubjectTag[2555] <- ds$SubjectTag[11]
  expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames=LoadDefaultOutcomeNames()), "The column 'SubjectTag' should not contain duplicated, but it does.")
})


test_that("outcomeNames -Missing argument entirely ", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)
  expect_error(ValidateOutcomeDataset(dsOutcome=ds), "The parameter for 'outcomeNames' should be passed, but was not.")
})


test_that("outcomeNames -Missing column", {
  ds <- LoadOutcomeFile()
  ds$SubjectTag <- CreateSubjectTag(ds$SubjectID, ds$Generation)
  expect_error(ValidateOutcomeDataset(dsOutcome=ds, outcomeNames="BAD"))#, "The parameter for 'outcomeNames' should be passed, but was not.")
})
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
