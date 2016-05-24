# Test for loss function analysis function
# 
# Author: Emilio L. Cano
###############################################################################


context("Loss Function Analysis")

lfaNominal <- ss.lfa(ss.data.bolts, "diameter", 
  lfa.Delta = 0.5, lfa.L0 = 0.001, lfa.Y0 = 10, lfa.output = "text")
lfaSmaller <- ss.lfa(ss.data.bolts, "diameter", 
  lfa.Delta = 0.5, lfa.L0 = 0.001, lfa.Y0 = 0, lfa.output = "text")
lfaLarger <- ss.lfa(ss.data.bolts, "diameter", 
  lfa.Delta = 0.5, lfa.L0 = 0.001, lfa.Y0 = Inf, lfa.output = "text")

testthat::test_that("Constant k is correctly computed",{
    expect_that(lfaNominal$lfa.k, equals(0.001/(0.5^2)))
    expect_that(lfaSmaller$lfa.k, equals(0.001/(0.5^2)))
    expect_that(lfaLarger$lfa.k, equals(0.001*(0.5^2)))
  })

