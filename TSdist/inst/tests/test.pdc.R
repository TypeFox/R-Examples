context("PDC ")

test_that("The PDC distance is calculated correctly", {
  
  # We depart from an example in the documentation of the pdc package to
  # check that the wrapper works correctly: 
  
  x1 <- sin(1:500)+rnorm(500,0,.1)
  x2 <- sin(1:500)+rnorm(500,0,.2)
  x3 <- sin(1:500)+rnorm(500,0,.3)
  x4 <- sin(1:500)+rnorm(500,0,.4)
  x5 <- rnorm(500,0,1)
  x6 <- rnorm(500,0,1)
  X <- cbind(x1, x2, x3, x4, x5, x6)  
  
  
  d <- pdcDist(X, 3)
  d <- as.matrix(d)
  
  expect_equal(PDCDistance(x1, x2, 3), d[2, 1])
  expect_equal(PDCDistance(x1, x3, 3), d[3, 1])
  expect_equal(PDCDistance(x1, x4, 3), d[4, 1])
  expect_equal(PDCDistance(x1, x5, 3), d[5, 1])
  expect_equal(PDCDistance(x1, x6, 3), d[6, 1])
  expect_equal(PDCDistance(x2, x3, 3), d[3, 2])
  expect_equal(PDCDistance(x2, x4, 3), d[4, 2])
  expect_equal(PDCDistance(x2, x5, 3), d[5, 2])
  expect_equal(PDCDistance(x2, x6, 3), d[6, 2])
  expect_equal(PDCDistance(x3, x4, 3), d[4, 3])
  expect_equal(PDCDistance(x3, x5, 3), d[5, 3])
  expect_equal(PDCDistance(x3, x6, 3), d[6, 3])
  expect_equal(PDCDistance(x4, x5, 3), d[5, 4])
  expect_equal(PDCDistance(x4, x6, 3), d[6, 4])
  expect_equal(PDCDistance(x5, x6, 3), d[6, 5])
  
  # If both m and t are null
  expect_equal(PDCDistance(x1, x2), as.numeric(pdcDist(cbind(x1, x2))))
  
  # If the original function outputs errors, the wrapper will return NA
  x1[1] <- "a"
  expect_equal(PDCDistance(x1, x2), NA)
  
  

})
