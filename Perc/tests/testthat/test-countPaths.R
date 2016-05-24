# to do 
# testthat if conf is not in the right format, raise an error. 
# testthat if maxlength is of right format  
# testthat output is of right format 
# testthat output is correct. 

context("Testing countPaths Function")

# test data generation
set.seed(1)
edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE), 
                        col2 = sample(letters[1:15], 200, replace = TRUE), 
                        stringsAsFactors = FALSE)
edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
testMatrix2 <- as.conflictmat(edgelist1)




test_that("return error for incorrect maxLength", {
  
  expect_error(countPaths(testMatrix2, maxLength = 1), "'maxLength' should be no smaller than 2.")
  expect_error(countPaths(testMatrix2, maxLength = 7), "'maxLength' should be no greater than 6.")
  expect_error(countPaths(testMatrix2, maxLength = 2.3), "'maxLength' needs to be an integer.")
})



