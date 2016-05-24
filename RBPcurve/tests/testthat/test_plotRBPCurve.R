context("plotRBPCurve")

test_that("plotRBPCurve", {
  for (i in tf) {
    plotRBPCurve(obj, cond.axis = i)
    plotRBPCurve(obj, add = TRUE, col = 2) 
  }
  dev.off()
  expect_error(plotRBPCurve(obj, add=TRUE), "plot.new has not been called yet")
})