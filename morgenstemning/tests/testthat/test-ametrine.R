context("ametrine")

test_that("output colors are equal to those produced by MATLAB code", {
  matlabOutput <- matrix(c(0.1176, 0.2353, 0.5882,
                           0.3137, 0.2745, 0.5948,
                           0.5098, 0.3137, 0.6013,
                           0.7059, 0.3529, 0.6078,
                           0.7712, 0.3464, 0.4902,
                           0.8366, 0.3399, 0.3725,
                           0.9020, 0.3333, 0.2549,
                           0.8889, 0.5098, 0.1699,
                           0.8758, 0.6863, 0.0850,
                           0.8627, 0.8627, 0.0000),
                         ncol=3, byrow=T)
  matlabCols <- rgb(matlabOutput)
  rCols <- ametrine(10)
  
  expect_equal(matlabCols, rCols)
})

test_that("invert option inverts correctly", {
  rcols <- ametrine(10)
  rcols.inv <- rev(rcols)
  rinvcols <- ametrine(10, invert=TRUE)
  
  expect_equal(rcols.inv, rinvcols)
})

test_that("minColor and maxColor arguments are handled correctly", {
  col.min <- "#FF0000"
  col.max <- "#00FF00"
  rcols <- ametrine(15, mincolor=col.min, maxcolor=col.max)
  
  expect_equal(rcols[1], col.min)
  expect_equal(rcols[length(rcols)], col.max)
})

test_that("alpha values are set correctly", {
  rcols <- ametrine(10, alpha=0.8)
  expect_equal(substr(rcols[2], 8, 9), "CC")
})
