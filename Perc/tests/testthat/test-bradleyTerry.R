context("Testing bradleyTerry function")


# test data generation
set.seed(1)
edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE), 
                        col2 = sample(letters[1:15], 200, replace = TRUE), 
                        stringsAsFactors = FALSE)
edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
testMatrix3 <- as.conflictmat(edgelist1)

test_that("Conflict Matrix does not meet Bradley-Terry assumption.", {
  testMatrix2 <- as.conflictmat(edgelist1)
  testMatrix2[1,] <- 0
  expect_error(bradleyTerry(testMatrix2),
               "Conflict Matrix does not meet Bradley-Terry assumption.  MLE does not exist.")

  })

test_that("output is a list of length 2", {
  bradleyTerryOutput <- bradleyTerry(testMatrix3)
  
  expect_output(str(bradleyTerryOutput), "List of 3")
})


test_that("outputs are correct", {
  
  expect_equal_to_reference(bradleyTerry(testMatrix3), file = "bradleyTerryoutput1.rds")
  
})