library(noncompliance)
context("Finding p-values")

test_that("Finding the p-values", {
  expect_equal(object = Get_pvalues_CONT(16, 1, 5, 1, 2, 8),
               expected = 0.002253486,
               tolerance = 1e-6)
  expect_equal(object = Get_pvalues_CONT(16, 1, 5, 1, 2, 8,
                                         useGLR = TRUE, justexactp = FALSE),
               expected = 0.005454856,
               tolerance = 1e-6)
})
