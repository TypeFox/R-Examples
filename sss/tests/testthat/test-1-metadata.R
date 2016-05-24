# Add comment
# 
# Author: Andrie
#----------------------------------------------------------------------------------

if(interactive()){
  library(testthat)
  sampleRoot <- "sss/tests/testthat/samples/sample-0"
} else {
  sampleRoot <- "samples/sample-0"
}
filenameSSS <- file.path(sampleRoot, "sample.sss")
filenameASC <- file.path(sampleRoot, "sample.asc")

context("Test metadata functions")

test_that("readSSSmetadata works", {
      
      test <- readSSSmetadata(filenameSSS)
      expect_is(test, "XMLDocumentContent")
      
    })

test_that("parseSSSmetadata works", {
      
      test <- parseSSSmetadata(readSSSmetadata(filenameSSS))
      rest <- structure(list(
              variables = structure(list(
                      ident = c("1", "2", "3", "4", "5", "6", "7", "99"), 
                      type = c("single", "multiple", "character", "multiple", 
                          "quantity", "logical", "single", "quantity"),
                      name = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q99"), 
                      label = c("Number of visits", "Attractions visited", "Other attractions visited", 
                          "Two favourite attractions visited", "Miles travelled", "Would come again", 
                          "When is that most likely to be", "Case weight"), 
                      positionStart = c(1, 2, 11, 41, 43, 46, 47, 48), 
                      positionFinish = c(1, 10, 40, 42, 45, 46, 47, 54), 
                      subfields = c("0", "0", "0", "2", "0", "0", "0", "0"), 
                      width = c(0, 0, 0, 1, 0, 0, 0, 0), 
                      hasValues = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)), 
                  .Names = c("ident", "type", "name", "label", "positionStart", "positionFinish", 
                      "subfields", "width", "hasValues"), 
                  row.names = c(NA, 8L), 
                  class = "data.frame"), 
              codes = structure(list(
                      ident = c("1", "1", "1", "2", "2", "2", "2", "2", "2", "3", 
                          "4", "4", "4", "4", "4", "4", "5", "5", "6", "7", "7", "7", "99"), 
                      code = c("1", "2", "3", "1", "2", "3", "4", "5", "9", NA, "1", "2", "3", "4", 
                          "5", "9", "500", "999", NA, "1", "2", "3", NA), 
                      codevalues = c("First visit", "Visited before within the year", "Visited before that", "Sherwood Forest", 
                          "Nottingham Castle", "\"Friar Tuck\" Restaurant", "\"Maid Marion\" Cafe", 
                          "Mining museum", "Other", NA, "Sherwood Forest", "Nottingham Castle", 
                          "\"Friar Tuck\" Restaurant", "\"Maid Marion\" Cafe", "Mining museum", 
                          "Other", "500 or more", "Not stated", NA, "Within 3 months", 
                          "Between 3 months and 1 year", "More than 1 years time", NA)), 
                  .Names = c("ident", "code", "codevalues"), 
                  row.names = c(NA, 23L), 
                  class = "data.frame")), 
              .Names = c("variables", "codes"))
      
      expect_is(test, "list")
      expect_equal(test$variables, rest$variables)
      expect_equal(test$codes, rest$codes)
      
    })

