library(testthat)
library(rpf)

context("ordinal.gamma")

# Example data from Agresti (1990, p. 21)
jobsat <- matrix(c(20,22,13,7,24,38,28,18,80,104,81,54,82,125,113,92), nrow=4, ncol=4)
expect_equal(ordinal.gamma(jobsat), .1265, tolerance=.001)
