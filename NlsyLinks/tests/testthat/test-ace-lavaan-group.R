options(digits=20)

###########
context("Lavaan")
###########
test_that("AceLavaanGroup -MathStandardized", {
  dsFull <- Links79PairExpanded #Start with the built-in data.frame in NlsyLinks
  oName_S1 <- "MathStandardized_S1" #Stands for Manifest1
  oName_S2 <- "MathStandardized_S2" #Stands for Manifest2
  dsFull <- dsFull[dsFull$RelationshipPath=='Gen2Siblings', ]
  
  dsGroupSummary <- RGroupSummary(dsFull, oName_S1, oName_S2)
#   rLevels <- dsGroupSummary[dsGroupSummary$Included, "R"]
  dsClean <- CleanSemAceDataset(dsDirty=dsFull, dsGroupSummary=dsGroupSummary, oName_S1=oName_S1, oName_S2=oName_S2)
#   dsClean
  
  ace <- AceLavaanGroup(dsClean)
  
  expectedASquared <- 0.621097662459077
  expectedCSquared  <- 0.210035603916836
  expectedESquared <- 0.168866733624087
  expectedCaseCount <- 8338 #8292
    
  expect_equal(object=ace@ASquared, expected=expectedASquared, scale=1)
  expect_equal(object=ace@CSquared, expected=expectedCSquared, scale=1)
  expect_equal(object=ace@ESquared, expected=expectedESquared, scale=1)
  expect_equal(object=ace@CaseCount, expected=expectedCaseCount, scale=1)
  expect_equal(object=slot(ace, "CaseCount"), expected=expectedCaseCount, scale=1)
  expect_true(object=slot(ace, "Unity"))
  expect_true(object=slot(ace, "WithinBounds"))
})
test_that("AceLavaanGroup -HeightZGenderAge", {
  dsFull <- Links79PairExpanded #Start with the built-in data.frame in NlsyLinks
  oName_S1 <- "HeightZGenderAge_S1"
  oName_S2 <- "HeightZGenderAge_S2"
  dsFull <- dsFull[dsFull$RelationshipPath=='Gen2Siblings', ]
  
  dsGroupSummary <- RGroupSummary(dsFull, oName_S1, oName_S2)
  #   rLevels <- dsGroupSummary[dsGroupSummary$Included, "R"]
  dsClean <- CleanSemAceDataset(dsDirty=dsFull, dsGroupSummary=dsGroupSummary, oName_S1=oName_S1, oName_S2=oName_S2)
  
  ace <- AceLavaanGroup(dsClean)
  
  expectedASquared <- 0.785694103632885
  expectedCSquared  <- 0.0320796445782708
  expectedESquared <- 0.182226251788844
  expectedCaseCount <- 5884
  
  expect_equal(object=ace@ASquared, expected=expectedASquared, scale=1)
  expect_equal(object=ace@CSquared, expected=expectedCSquared, scale=1)
  expect_equal(object=ace@ESquared, expected=expectedESquared, scale=1)
  expect_equal(object=ace@CaseCount, expected=expectedCaseCount, scale=1)
  expect_true(object=slot(ace, "Unity"))
  expect_true(object=slot(ace, "WithinBounds"))
})
# test_that("AceLavaanGroup -WeightStandardizedForAge19To25", {
#   dsFull <- Links79PairExpanded #Start with the built-in data.frame in NlsyLinks
#   oName_S1 <- "WeightStandardizedForAge19To25_S1"
#   oName_S2 <- "WeightStandardizedForAge19To25_S2"
#   dsFull <- dsFull[dsFull$RelationshipPath=='Gen2Siblings', ]
#   
#   dsGroupSummary <- RGroupSummary(dsFull, oName_S1, oName_S2)
# #   rLevels <- dsGroupSummary[dsGroupSummary$Included, "R"]
#   dsClean <- CleanSemAceDataset(dsDirty=dsFull, dsGroupSummary=dsGroupSummary, oName_S1=oName_S1, oName_S2=oName_S2)
#   
#   ace <- AceLavaanGroup(dsClean)
#   
#   expectedASquared <- 0.68859786571166337 #.68764112310158664876 #0.687801119966999
#   expectedCSquared  <- 2.0347745721123465e-17 #1.6026541226325013592e-17 #7.10884537223939e-15
#   expectedESquared <- 0.31140213428833668 #.31235887689841324022 #0.312198880032994
#   expectedCaseCount <- 3479 #3478
#   
#   expect_equal(object=ace@ASquared, expected=expectedASquared, scale=1)
#   expect_equal(object=ace@CSquared, expected=expectedCSquared, scale=1)
#   expect_equal(object=ace@ESquared, expected=expectedESquared, scale=1)
#   expect_equal(object=ace@CaseCount, expected=expectedCaseCount, scale=1)
#   expect_true(object=slot(ace, "Unity"))
#   expect_true(object=slot(ace, "WithinBounds"))
# })

# stringr::str_c(ace@ASquared)
# stringr::str_c(ace@CSquared)
# stringr::str_c(ace@ESquared)
# stringr::str_c(ace@CaseCount)
