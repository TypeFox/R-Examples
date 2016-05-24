# to do 
# testthat if conf is not in the right format, raise an error. -- checked
# testthat if maxlength is of right format  -- checked
# testthat output is of right format -- checked
# testthat output is correct. -- checked

context("Testing findAllPaths Function")

# test data generation
set.seed(1)
edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE), 
                        col2 = sample(letters[1:15], 200, replace = TRUE), 
                        stringsAsFactors = FALSE)
edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
testMatrix2 <- as.conflictmat(edgelist1)


test_that("return error for incorrect maxLength", {
  
  expect_error(findAllPaths(testMatrix2, maxLength = 1), "'maxLength' should be no smaller than 2.")
  expect_error(findAllPaths(testMatrix2, maxLength = 7), "'maxLength' should be no greater than 6.")
  expect_error(findAllPaths(testMatrix2, maxLength = 2.3), "'maxLength' needs to be an integer.")
})

test_that("output is a list of the right structure", {
  outputlist_maxLength2 <- countPaths(testMatrix2, maxLength = 2)
  outputlist_maxLength4 <- countPaths(testMatrix2, maxLength = 4)
  
  expect_is(outputlist_maxLength2, "list")
  expect_is(outputlist_maxLength4, "list")
  
  expect_equal(length(outputlist_maxLength2), 1)
  expect_equal(length(outputlist_maxLength4), 3)
})


test_that("outputs are correct", {
  
  expect_equal_to_reference(countPaths(testMatrix2, maxLength = 2), file = "countPathsOutput1.rds")
  expect_equal_to_reference(countPaths(testMatrix2, maxLength = 3), file = "countPathsOutput2.rds")
  expect_equal_to_reference(countPaths(testMatrix2, maxLength = 4), file = "countPathsOutput3.rds")
})

