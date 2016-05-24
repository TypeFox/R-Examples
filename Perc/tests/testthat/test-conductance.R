# to do 
# testthat if conf is not in the right format, raise an error. -- checked
# testthat if maxLength is of right format -- checked
# testthat output is a list of two elements. -- checked
# testthat output is correct. -- checked


context("Testing Conductance Function")

# test data generation
set.seed(1)
edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE), 
                        col2 = sample(letters[1:15], 200, replace = TRUE), 
                        stringsAsFactors = FALSE)
edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
testMatrix2 <- as.conflictmat(edgelist1)


test_that("return error for incorrect maxLength", {
  
  expect_error(conductance(testMatrix2, maxLength = 1), "'maxLength' should be an integer greater than 1 and less than 7.")
  expect_error(conductance(testMatrix2, maxLength = 2.3), "'maxLength' needs to be an integer.")
  
})

test_that("output is a list of length 2", {
  conductanceOutput <- conductance(testMatrix2, maxLength = 2, strict = FALSE)
  expect_is(conductanceOutput, "list")
  expect_equal(length(conductanceOutput), 2)
  expect_output(str(conductanceOutput), "List of 2")
})



test_that("outputs are correct", {
 
  expect_equal_to_reference(conductance(testMatrix2, maxLength = 2, strict = FALSE), file = "conductanceOutput1.rds")
  expect_equal_to_reference(conductance(testMatrix2, maxLength = 3, strict = FALSE), file = "conductanceOutput2.rds")
  expect_equal_to_reference(conductance(testMatrix2, maxLength = 4, strict = FALSE), file = "conductanceOutput3.rds")
})