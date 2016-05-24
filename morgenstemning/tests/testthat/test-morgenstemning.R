context("morgenstemning")

test_that("output colors are equal to those produced by MATLAB code", {
  matlabOutput <- matrix(c(0, 0, 0,
                           0.0183, 0.1046, 0.1662,
                           0.0759, 0.1922, 0.3338,
                           0.3142, 0.1919, 0.4467,
                           0.6283, 0.1301, 0.4497,
                           0.8788, 0.2421, 0.3711,
                           0.9683, 0.6527, 0.1342,
                           0.9896, 0.9388, 0.1076,
                           0.9982, 0.9978, 0.6244,
                           1.0000, 1.0000, 1.0000),
                    ncol=3, byrow=T)
  matlabCols <- rgb(matlabOutput)
  rCols <- morgenstemning(10)
  
  expect_equal(matlabCols, rCols)
})

test_that("invert option inverts correctly", {
  rcols <- morgenstemning(10)
  rcols.inv <- rev(rcols)
  rinvcols <- morgenstemning(10, invert=TRUE)
  
  expect_equal(rcols.inv, rinvcols)
})

test_that("minColor and maxColor arguments are handled correctly", {
  col.min <- "#FF0000"
  col.max <- "#00FF00"
  rcols <- morgenstemning(15, mincolor=col.min, maxcolor=col.max)
  
  expect_equal(rcols[1], col.min)
  expect_equal(rcols[length(rcols)], col.max)
})

test_that("alpha values are set correctly", {
  rcols <- morgenstemning(10, alpha=0.8)
  expect_equal(substr(rcols[2], 8, 9), "CC")
})
