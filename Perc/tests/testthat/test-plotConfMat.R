# to do 
# testthat if conf is not in the right format, raise an error. -- checked
# testthat output is a list of 45 elements.                    -- checked
# testthat output is correct.                                  -- checked


context("Testing plotConfmat Function")

# test data generation
set.seed(1)
edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE), 
                        col2 = sample(letters[1:15], 200, replace = TRUE), 
                        stringsAsFactors = FALSE)
edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
testMatrix2 <- as.conflictmat(edgelist1)

# tests
test_that("input 'conf' is a matrix.", {
  expect_error(plotConfmat(edgelist1),
               "conf.mat should be a matrix.")
})

test_that("output is a list of length 45", {
  expect_output(str(plotConfmat(testMatrix2)), "List of 45")
  expect_output(
    str(plotConfmat(testMatrix2, ordering = TRUE, labels = TRUE)), 
    "List of 45"
    )
})

test_that("outputs are correct", {
  
  expect_equal_to_reference(
    plotConfmat(testMatrix2, ordering = TRUE, labels = TRUE), 
    file = "plotConfMatOutput1.rds"
    )

})