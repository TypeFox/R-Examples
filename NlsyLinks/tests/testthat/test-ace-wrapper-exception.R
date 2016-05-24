
###########
context("Ace Wrapper Exceptions")
###########
test_that("AceUnivariate -NULL method", {
  dsOutcomes <- ExtraOutcomes79
  dsOutcomes$SubjectTag <- CreateSubjectTag(subjectID=dsOutcomes$SubjectID,generation=dsOutcomes$Generation)
  dsDF <- CreatePairLinksDoubleEntered(outcomeDataset=dsOutcomes, linksPairDataset=Links79Pair, outcomeNames=c("MathStandardized", "WeightZGenderAge"))
  
  expect_error(
    AceUnivariate(method=NULL, dataSet=dsOutcomes, oName_S1="MathStandardized_S1", oName_S2="MathStandardized_S2"),
    "The method argument must contain exactly one element when calling the AceUnivariate function.  It contained 0 elements."
  )
})
test_that("AceUnivariate -blank method", {
  dsOutcomes <- ExtraOutcomes79
  dsOutcomes$SubjectTag <- CreateSubjectTag(subjectID=dsOutcomes$SubjectID,generation=dsOutcomes$Generation)
  dsDF <- CreatePairLinksDoubleEntered(outcomeDataset=dsOutcomes, linksPairDataset=Links79Pair, outcomeNames=c("MathStandardized", "WeightZGenderAge"))
  
  expect_error(
    AceUnivariate(method="", dataSet=dsOutcomes, oName_S1="MathStandardized_S1", oName_S2="MathStandardized_S2"),
    "The method argument must contain exactly one element when calling the AceUnivariate function.  It was blank."
  )
})
test_that("AceUnivariate -unrecognized name method", {
  dsOutcomes <- ExtraOutcomes79
  dsOutcomes$SubjectTag <- CreateSubjectTag(subjectID=dsOutcomes$SubjectID,generation=dsOutcomes$Generation)
  dsDF <- CreatePairLinksDoubleEntered(outcomeDataset=dsOutcomes, linksPairDataset=Links79Pair, outcomeNames=c("MathStandardized", "WeightZGenderAge"))
  
  expect_error(
    AceUnivariate(method="ddddd", dataSet=dsOutcomes, oName_S1="MathStandardized_S1", oName_S2="MathStandardized_S2"),
    "The method argument, 'ddddd' was not recognized as a valid option to the AceUnivariate function."
  )
})
test_that("AceUnivariate -multiple methods passed", {
  dsOutcomes <- ExtraOutcomes79
  dsOutcomes$SubjectTag <- CreateSubjectTag(subjectID=dsOutcomes$SubjectID,generation=dsOutcomes$Generation)
  dsDF <- CreatePairLinksDoubleEntered(outcomeDataset=dsOutcomes, linksPairDataset=Links79Pair, outcomeNames=c("MathStandardized", "WeightZGenderAge"))
  
  expect_error(
    AceUnivariate(method=c("DeFriesFulkerMethod1", "DeFriesFulkerMethod3"), dataSet=dsOutcomes, oName_S1="MathStandardized_S1", oName_S2="MathStandardized_S2"),
    "The method argument must contain exactly one element when calling the AceUnivariate function.  It contained 2 elements."
  )
})
