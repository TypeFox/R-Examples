#require(testthat); library("gapfill", lib.loc = "../../../lib/")
context("test-Predict")

test_that("Predict",{
  a <- ArrayAround(ndvi, 
                   mp = IndexOneFour(which(is.na(ndvi))[10], dim(ndvi)),
                   size = c(10,10,1,1))
  expect_equal(Predict(a), 0.6745)

  ## arg: nTargetImage
  expect_equal(Predict(a, nTargetImage = 203), NA)
  aa <- a
  aa[,,2,1] <- NA
  aa[1:2,1:2,2,1] <- a[1:2,1:2,2,1]
  expect_equal(Predict(aa), NA)

  ## arg: nImages
  expect_equal(Predict(a, nImages = 7), NA)
  aa <- a
  aa[,,1,1:2] <- NA
  aa[,,2,1] <- NA
  expect_equal(Predict(aa), NA)

  ## arg: nQuantile
  expect_equal(Predict(a, nQuant = 3), 0.6229)
  

  ## EstimateQuantile
  expect_equal(round(EstimateQuantile(a, attr(a, "mp"), 2), 4),
               0.8849)


})
