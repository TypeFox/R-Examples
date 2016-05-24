context("DTW")

test_that("The DTW function is calculated correctly", {
  
  # We use an example in the documentation of the dtw package to test that the wrapper works correctly: 
  
  idx <- seq(0, 6.28, len=100);
  query <- sin(idx) + runif(100) / 10;
  reference <- cos(idx)
  
  d <- dtw(query,reference)$d
  expect_equal(DTWDistance(query, reference), d)

  # If an error is thrown by the dtw function, DTWDistance will return NA.
  # For example, if there are missing values in the series:
  query[1] <- NA
  expect_equal(DTWDistance(query, reference), NA)
  
  # For specific errors of the dtw function access the documentation of the 
  # dtw package.
   
})
