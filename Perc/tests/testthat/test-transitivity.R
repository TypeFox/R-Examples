# to do 
# testthat if conf is not in the right format, raise an error.  -- checked
# testthat output is a list of four elements.                   -- checked
# testthat output is correct. 


context("Testing transitivity Function")

# test data generation
set.seed(1)
edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE), 
                        col2 = sample(letters[1:15], 200, replace = TRUE), 
                        stringsAsFactors = FALSE)
edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
testMatrix2 <- as.conflictmat(edgelist1)

# tests

test_that("output is a list of length 4", {
  expect_is(transitivity(testMatrix2, strict = FALSE), "list")
  expect_equal(length(transitivity(testMatrix2)), 4)
})

test_that("outputs are correct", {
  expect_equal_to_reference(transitivity(testMatrix2, strict = FALSE), file = "transitivityOutput1.rds")
})