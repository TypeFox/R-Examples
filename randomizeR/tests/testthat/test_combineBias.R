###################################################################
# ----------------------------------------------------------------#
# Tests for objects from the bias class and associated functions  #
# ----------------------------------------------------------------#
###################################################################
context("combineBias")

test_that("returns valid object", {
  method <- sample(c("exact", "sim"), 1)
  cB <- chronBias(type = sample(c("linT", "stepT", "logT"), 1), 
                         theta = runif(1), 
                         saltus = sample(10, 1),
                         method = method)
  sB <- selBias(type = sample(c("CS", "DS"), 1),
                eta = runif(1),
                method = method)
  expect_is(combineBias(sB, cB), "combinedBias")   
  expect_error(combineBias(sB, sB))
  expect_error(combineBias(cB, cB))
  expect_error(combineBias(cB, sB))
})    
    

