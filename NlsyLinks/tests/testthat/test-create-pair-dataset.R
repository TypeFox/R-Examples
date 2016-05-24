
#options(digits=17)
#expect_equal(object=10.01, expected=10, tolerance=.002, scale=1) #Absolute Difference
#expect_equal(object=10.01, expected=10, tolerance=.002, scale=NULL) #Relative Difference

LoadPairFile <- function( ) {
#   directory <- "F:/Projects/Nls/Links2011/Analysis/Df/2012-01-13"
#   pathLinks <- file.path(directory, "Links2011V28.csv")
  #dsLinks <- read.csv(pathLinks)
  return( Links79Pair )
}
LoadOutcomeFile <- function( ) {
#   directory <- "F:/Projects/Nls/Links2011/Analysis/Df/2012-01-13"
#   pathDv <-  file.path(directory, "BMI_Sex_Intell.csv")
#   dsDv <- read.csv(pathDv)
  return( ExtraOutcomes79 )
}
LoadDefaultOutcomeNames <- function( ) {
  return( c("HeightZGenderAge") )
}

###########
context("CreatePairLinksDoubleEntered")
###########
test_that("CreatePairLinksDoubleEntered -Normal Scenario", {
  dsLinks <- LoadPairFile()
  dsLinks <- dsLinks[dsLinks$RelationshipPath=='Gen2Siblings', ]
  
  dsOutcomes <- LoadOutcomeFile()
  dsOutcomes$SubjectTag <- CreateSubjectTag(subjectID=dsOutcomes$SubjectID, generation=dsOutcomes$Generation)
  dsLinksWithExtraOutcome <- CreatePairLinksDoubleEntered(outcomeNames=LoadDefaultOutcomeNames(), outcomeDataset=dsOutcomes, linksPairDataset=dsLinks)
  expect_equal(nrow(dsLinksWithExtraOutcome), 22176, info="The number of rows in the pairs links should be correct.")
  expect_equal(ncol(dsLinksWithExtraOutcome), 7, info="The number of columns in the pairs links should be correct.")  
  
  expectedColumnNames <- c("SubjectTag_S1", "SubjectTag_S2", "ExtendedID", "R", "RelationshipPath", "HeightZGenderAge_S1", "HeightZGenderAge_S2")
  actualColumnNames <- colnames(dsLinksWithExtraOutcome)
  expect_equal(actualColumnNames, expectedColumnNames, info="The column names, and their order, should be correct.")
  
  expect_equal(mean(dsLinksWithExtraOutcome$HeightZGenderAge_S1, na.rm=T), 0, tolerance=.05, scale=1)
  expect_equal(mean(dsLinksWithExtraOutcome$HeightZGenderAge_S1, na.rm=T), mean(dsLinksWithExtraOutcome$HeightZGenderAge_S2, na.rm=T))#, tolerance=.01, scale=1)
  
#   expect_equal(mean(dsLinksWithExtraOutcome$Weight_S1, na.rm=T), 161.949, tolerance=1e-4, scale=1)
#   expect_equal(mean(dsLinksWithExtraOutcome$Weight_S1, na.rm=T), mean(dsLinksWithExtraOutcome$Weight_S2, na.rm=T))#, tolerance=.01, scale=1)  

  #The following aren't tested against meaningful values, but they do provide some regression testing.
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$SubjectTag_S1), na.rm=T), 13164976703)
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$SubjectTag_S1), na.rm=T), sum(as.numeric(dsLinksWithExtraOutcome$SubjectTag_S2, na.rm=T)))  
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$ExtendedID), na.rm=T), 131575156)
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$R), na.rm=T), 9247.5) #9235.75)# 8813.5)  
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$RelationshipPath), na.rm=T), 44352)
})
test_that("CreatePairLinksDoubleEntered -Normal Scenario 2 sibs", {
  dsExpected <- data.frame(
    SubjectTag_S1=c(101, 102),
    SubjectTag_S2=c(102, 101), 
    ExtendedID=c(1, 1),
    R=c(.5, .5), 
    RelationshipPath=rep("Gen2Siblings", 2),
    DV1_S1=c(11, 12),
    DV2_S1=c(21, 22),
    DV1_S2=c(12, 11),
    DV2_S2=c(22, 21)
  )
  dsSingleLinks <- data.frame(ExtendedID=c(1), SubjectTag_S1=c(101), SubjectTag_S2=c(102), R=c(.5), RelationshipPath=rep("Gen2Siblings", 1))
  dsSingleOutcomes <- data.frame(SubjectTag=c(101, 102), DV1=c(11, 12), DV2=c(21, 22))
  dsDouble <- CreatePairLinksDoubleEntered(outcomeDataset=dsSingleOutcomes, linksPairDataset=dsSingleLinks, outcomeNames=c("DV1", "DV2"), validateOutcomeDataset=T)
  
  expect_equal(object=dsDouble, expected=dsExpected)
})
test_that("CreatePairLinksDoubleEntered -Normal Scenario 3 sibs", {
  dsExpected <- data.frame(
    SubjectTag_S1=c(101, 101, 102, 102, 103, 103), 
    SubjectTag_S2=c(102, 103, 103, 101, 101, 102), 
    ExtendedID=rep(1, 6), 
    R=c(.5, .25, .25, .5, .25, .25), 
    RelationshipPath=rep("Gen2Siblings", 6),
    DV1_S1=c(11, 11, 12, 12, 13, 13),
    DV2_S1=c(21, 21, 22, 22, 23, 23),    
    DV1_S2=c(12, 13, 13, 11, 11, 12),    
    DV2_S2=c(22, 23, 23, 21, 21, 22)    
  )
  
  dsSingleLinks <- data.frame(ExtendedID=c(1, 1, 1), SubjectTag_S1=c(101, 101, 102), SubjectTag_S2=c(102, 103, 103), R=c(.5, .25, .25), RelationshipPath=rep("Gen2Siblings", 3))
  dsSingleOutcomes <- data.frame(SubjectTag=c(101, 102, 103), DV1=c(11, 12, 13), DV2=c(21, 22, 23))
  dsDouble <- CreatePairLinksDoubleEntered(outcomeDataset=dsSingleOutcomes, linksPairDataset=dsSingleLinks, outcomeNames=c("DV1", "DV2"), validateOutcomeDataset=T)
  
  expect_equal(object=dsDouble, expected=dsExpected)
})
test_that("CreatePairLinksDoubleEntered -Normal Scenario 2 familes", {
  #dsExpected <- data.frame(SubjectTag_S1=c(101, 101, 102, 201, 102, 103, 103, 202), SubjectTag_S2=c(102, 103, 103, 202, 101, 101, 102, 201), ExtendedID=c(1,1,1,2, 1,1,1,2), R=c(.5, .25, .25, .5, .5, .25, .25, .5), RelationshipPath=rep("Gen2Siblings", 8))
  dsExpected <- data.frame(
    SubjectTag_S1=c(101, 101, 102, 201, 102, 103, 103, 202), 
    SubjectTag_S2=c(102, 103, 103, 202, 101, 101, 102, 201), 
    ExtendedID=c(1,1,1,2, 1,1,1,2), 
    R=c(.5, .25, .25, .5, .5, .25, .25, .5), 
    RelationshipPath=rep("Gen2Siblings", 8),
    DV1_S1=c(11, 11, 12, 41, 12, 13, 13, 42),
    DV2_S1=c(21, 21, 22, 51, 22, 23, 23, 52),    
    DV1_S2=c(12, 13, 13, 42, 11, 11, 12, 41),    
    DV2_S2=c(22, 23, 23, 52, 21, 21, 22, 51)    
  )  
 
  dsSingleLinks <- data.frame(ExtendedID=c(1, 1, 1, 2), SubjectTag_S1=c(101, 101, 102, 201), SubjectTag_S2=c(102, 103, 103, 202), R=c(.5, .25, .25, .5), RelationshipPath=rep("Gen2Siblings", 4))
  dsSingleOutcomes <- data.frame(SubjectTag=c(101, 102, 103, 201, 202), DV1=c(11, 12, 13, 41, 42), DV2=c(21, 22, 23, 51, 52))
  dsDouble <- CreatePairLinksDoubleEntered(outcomeDataset=dsSingleOutcomes, linksPairDataset=dsSingleLinks, outcomeNames=c("DV1", "DV2"), validateOutcomeDataset=T)
 
  expect_equal(object=dsDouble, expected=dsExpected)
})


###########
context("CreatePairLinksSingleEntered")
###########
test_that("CreatePairLinksSingleEntered -Normal Scenario", {
  dsLinks <- LoadPairFile()
  dsLinks <- dsLinks[dsLinks$RelationshipPath=='Gen2Siblings', ]
  dsOutcomes <- LoadOutcomeFile()
  dsOutcomes$SubjectTag <- CreateSubjectTag(subjectID=dsOutcomes$SubjectID, generation=dsOutcomes$Generation)
  dsLinksWithExtraOutcome <- CreatePairLinksSingleEntered(outcomeNames=LoadDefaultOutcomeNames(), outcomeDataset=dsOutcomes, linksPairDataset=dsLinks)
  expect_equal(nrow(dsLinksWithExtraOutcome), 11088, info="The number of rows in the pairs links should be correct.")
  expect_equal(ncol(dsLinksWithExtraOutcome), 7, info="The number of columns in the pairs links should be correct.")  
  
  expectedColumnNames <- c("SubjectTag_S1", "SubjectTag_S2", "ExtendedID", "R", "RelationshipPath", "HeightZGenderAge_S1", "HeightZGenderAge_S2")
  actualColumnNames <- colnames(dsLinksWithExtraOutcome)
  expect_equal(actualColumnNames, expectedColumnNames, info="The column names, and their order, should be correct.")
  
  expect_equal(mean(dsLinksWithExtraOutcome$HeightZGenderAge_S1, na.rm=T), 0, tolerance=.04, scale=1)  
  expect_equal(mean(dsLinksWithExtraOutcome$HeightZGenderAge_S2, na.rm=T), 0, tolerance=.06, scale=1)
    
#   expect_equal(mean(dsLinksWithExtraOutcome$Weight_S1, na.rm=T),  167.03616780045351, tolerance=1e-15, scale=1)
#   expect_equal(mean(dsLinksWithExtraOutcome$Weight_S2, na.rm=T), 155.87657328461225, tolerance=1e-15, scale=1)
  
  #The following aren't tested against meaningful values, but they do provide some regression testing.
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$SubjectTag_S1), na.rm=T), 6582479134)
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$SubjectTag_S2), na.rm=T), 6582497569)  
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$ExtendedID), na.rm=T), 65787578)
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$R), na.rm=T), 4623.75)  
  expect_equal(sum(as.numeric(dsLinksWithExtraOutcome$RelationshipPath), na.rm=T), 22176)
})

###########
context("CreatePairLinksDoubleEnteredWithNoOutcomes")
###########
test_that("CreatePairLinksDoubleEnteredWithNoOutcomes -Normal Scenario 2 sibs", {
  dsExpected <- data.frame(SubjectTag_S1=c(101, 102), SubjectTag_S2=c(102, 101), ExtendedID=c(1, 1),R=c(.5, .5), RelationshipPath=rep("Gen2Siblings", 2))
  dsSingle <- data.frame(ExtendedID=c(1), SubjectTag_S1=c(101), SubjectTag_S2=c(102), R=c(.5), RelationshipPath=rep("Gen2Siblings", 1))
  dsDouble <- CreatePairLinksDoubleEnteredWithNoOutcomes(linksPairDataset=dsSingle)
  expect_equal(object=dsDouble, expected=dsExpected)
  #linksPairDataset <- dsSingle
})
test_that("CreatePairLinksDoubleEnteredWithNoOutcomes -Normal Scenario 3 sibs", {
  dsExpected <- data.frame(SubjectTag_S1=c(101, 101, 102, 102, 103, 103), SubjectTag_S2=c(102, 103, 103, 101, 101, 102), ExtendedID=rep(1, 6), R=c(.5, .25, .25, .5, .25, .25), RelationshipPath=rep("Gen2Siblings", 6))
  dsSingle <- data.frame(ExtendedID=c(1, 1, 1), SubjectTag_S1=c(101, 101, 102), SubjectTag_S2=c(102, 103, 103), R=c(.5, .25, .25), RelationshipPath=rep("Gen2Siblings", 3))
  dsDouble <- CreatePairLinksDoubleEnteredWithNoOutcomes(linksPairDataset=dsSingle)
  expect_equal(object=dsDouble, expected=dsExpected)
})
test_that("CreatePairLinksDoubleEnteredWithNoOutcomes -Normal Scenario 2 familes", {
  dsExpected <- data.frame(SubjectTag_S1=c(101, 101, 102, 201, 102, 103, 103, 202), SubjectTag_S2=c(102, 103, 103, 202, 101, 101, 102, 201), ExtendedID=c(1,1,1,2, 1,1,1,2), R=c(.5, .25, .25, .5, .5, .25, .25, .5), RelationshipPath=rep("Gen2Siblings", 8))
  dsSingle <- data.frame(ExtendedID=c(1, 1, 1, 2), SubjectTag_S1=c(101, 101, 102, 201), SubjectTag_S2=c(102, 103, 103, 202), R=c(.5, .25, .25, .5), RelationshipPath=rep("Gen2Siblings", 4))
  dsDouble <- CreatePairLinksDoubleEnteredWithNoOutcomes(linksPairDataset=dsSingle)
  expect_equal(object=dsDouble, expected=dsExpected)
})


###########
context("ValidatePairLinks")
###########
test_that("CreatePairLinksDoubleEnteredWithNoOutcomes -Short Scenario", {
  dsLinks <- LoadPairFile()
  expect_true(ValidatePairLinks(dsLinks))
})

test_that("Zero rows", {
  dsLinks <- LoadPairFile()
  dsLinks <- dsLinks[0,]
  expect_error(ValidatePairLinks(dsLinks), "The linksPair data frame should have at least one row, but does not.")
})

test_that("Bad SubjectTag_S1", {
  dsLinks <- LoadPairFile()
  expect_true(ValidatePairLinks(dsLinks))
  colnames(dsLinks)[colnames(dsLinks)=="SubjectTag_S1"] <- "Bad"
  expect_error(ValidatePairLinks(dsLinks), "The column 'SubjectTag_S1' should exist in the linksPair file, but does not.")
})

test_that("Bad SubjectTag_S2", {
  dsLinks <- LoadPairFile()
  expect_true(ValidatePairLinks(dsLinks))
  colnames(dsLinks)[colnames(dsLinks)=="SubjectTag_S2"] <- "Bad"
  expect_error(ValidatePairLinks(dsLinks), "The column 'SubjectTag_S2' should exist in the linksPair file, but does not.")
})

test_that("Bad R", {
  dsLinks <- LoadPairFile()
  expect_true(ValidatePairLinks(dsLinks))
  colnames(dsLinks)[colnames(dsLinks)=="R"] <- "Bad"
  expect_error(ValidatePairLinks(dsLinks), "The column 'R' should exist in the linksPair file, but does not.")
})

###########
context("ValidatePairLinksAreSymmetric")
###########
test_that("ValidatePairLinksAreSymmetric -Normal Scenario", {
  dsLinks <- LoadPairFile()
  dsDouble <- CreatePairLinksDoubleEnteredWithNoOutcomes(linksPairDataset=dsLinks)
  expect_true(ValidatePairLinksAreSymmetric(dsDouble))
})
test_that("ValidatePairLinksAreSymmetric -Short Scenario", {
  dsDouble <- data.frame(
    SubjectTag_S1=c(101, 101, 102, 102, 103, 103), 
    SubjectTag_S2=c(102, 103, 103, 101, 101, 102), 
    ExtendedID=rep(1, 6), 
    R=c(.5, .25, .25, .5, .25, .25), 
    RelationshipPath=rep("Gen2Siblings", 6),
    DV1_S1=c(11, 11, 12, 12, 13, 13),
    DV2_S1=c(21, 21, 22, 22, 23, 23),    
    DV1_S2=c(12, 13, 13, 11, 11, 12),    
    DV2_S2=c(22, 23, 23, 21, 21, 22)    
    )
  expect_true(ValidatePairLinksAreSymmetric(dsDouble))
})

test_that("ValidatePairLinksAreSymmetric -not doubled 1", {
  dsLinks <- LoadPairFile() 
  expect_error(ValidatePairLinksAreSymmetric(dsLinks), label="The 'linksPair' dataset doesn't appear to be double-entered & symmetric.  The reciprocal of (SubjectTag_S1, SubjectTag_S2, R)=(101, 102, 0.5) was found 0 time(s).")
})
test_that("ValidatePairLinksAreSymmetric -not doubled 2", {
  dsLinks <- data.frame(ExtendedID=c(1), SubjectTag_S1=c(101), SubjectTag_S2=c(102), R=c(.5), RelationshipPath=rep("Gen2Siblings", 2))
  expect_error(ValidatePairLinksAreSymmetric(dsLinks), label="The 'linksPair' dataset doesn't appear to be double-entered & symmetric.  The reciprocal of (SubjectTag_S1, SubjectTag_S2, R)=(101, 102, 0.5) was found 0 time(s).")
})
test_that("ValidatePairLinksAreSymmetric -Assymetric Scenario sho", {
  dsLinks <- LoadPairFile()
  dsLinks$R[2] <- .87
  expect_error(ValidatePairLinksAreSymmetric(dsLinks), label="The 'linksPair' dataset doesn't appear to be double-entered & symmetric.  The reciprocal of (SubjectTag_S1, SubjectTag_S2, R)=(201, 202, 0.5) was found 0 time(s).")
})
