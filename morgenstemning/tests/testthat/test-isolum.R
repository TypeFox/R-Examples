context("isolum")

test_that("output colors are equal to those produced by MATLAB code", {
  matlabOutput <- matrix(c(0.3184, 0.6723, 0.8669,
                           0.4639, 0.6369, 0.8171,
                           0.5793, 0.5998, 0.7648,
                           0.6960, 0.5871, 0.6696,
                           0.8013, 0.5785, 0.5541,
                           0.8996, 0.5439, 0.4413,
                           0.9921, 0.4859, 0.3239,
                           0.9366, 0.6123, 0.2586,
                           0.8784, 0.7205, 0.1760,
                           0.8170, 0.8170, 0.0000),
                    ncol=3, byrow=T)
  matlabCols <- rgb(matlabOutput)
  rCols <- isolum(10)
  
  expect_equal(matlabCols, rCols)
})

test_that("invert option inverts correctly", {
  rcols <- isolum(10)
  rcols.inv <- rev(rcols)
  rinvcols <- isolum(10, invert=TRUE)
  
  expect_equal(rcols.inv, rinvcols)
})

test_that("minColor and maxColor arguments are handled correctly", {
  col.min <- "#FF0000"
  col.max <- "#00FF00"
  rcols <- isolum(15, mincolor=col.min, maxcolor=col.max)
  
  expect_equal(rcols[1], col.min)
  expect_equal(rcols[length(rcols)], col.max)
})

test_that("alpha values are set correctly", {
  rcols <- isolum(10, alpha=0.8)
  expect_equal(substr(rcols[2], 8, 9), "CC")
})
