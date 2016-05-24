# to do 
# testthat if matrix is not in the right format, raise an error. -- checked
# testthat output is a dataframe. 
# testthat output is correct.

context("Testing Converters")

set.seed(1)
edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE), 
                        col2 = sample(letters[1:15], 200, replace = TRUE), 
                        stringsAsFactors = FALSE)
edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
sampleMatrix <- as.conflictmat(edgelist1)
testMatrix <- conductance(sampleMatrix, 2)[[2]]

test_that("input 'conf' is of 'conf.mat'", {
  
  expect_error(valueConverter(edgelist1),
               "Only matrix is accepted as input.")
  expect_error(dyadicLongConverter(edgelist1),
               "Only matrix is accepted as input.")
  expect_error(individualDomProb(edgelist1),
               "Only matrix is accepted as input.")
})


test_that("outputs are of correct format for valueConverter", {
  vcOutput <- valueConverter(testMatrix)
  expect_is(vcOutput, "matrix")
  expect_equal(dim(vcOutput), dim(testMatrix))
})


test_that("outputs are of correct format for dyadicLongConverter", {
  longData <- dyadicLongConverter(testMatrix)
  expect_output(str(longData), "data.frame")
  expect_equal(ncol(longData), 5)
})

test_that("outputs are of correct format for individualDomProb", {
  iwOutput <- individualDomProb(testMatrix)
  expect_output(str(iwOutput), "data.frame")
  expect_equal(ncol(iwOutput), 3)
})

test_that("outputs are correct", {
  
  expect_equal_to_reference(valueConverter(testMatrix), file = "valueConverterOutput.rds")
  expect_equal_to_reference(dyadicLongConverter(testMatrix), file = "dyadicLongConverterOutput.rds")
  expect_equal_to_reference(individualDomProb(testMatrix), file = "individualDomProbOutput.rds")
})
