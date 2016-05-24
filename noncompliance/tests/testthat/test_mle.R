library(noncompliance)
context("Finding MLEs")

test_that("Finding MLE under sharp H0 for Compliers", {
  expect_equal(object = FindMLE_CONT_H0_hypergeoR(7, 4, 2, 3, 1, 5)$qH0,
               expected = 0.02000476,
               tolerance = 1e-6)
  expect_equal(object = FindMLE_CONT_H0_hypergeoR(16, 1, 5, 1, 2, 8)$qH0,
               expected = 0.0007159563,
               tolerance = 1e-6)
  expect_equal(object = FindMLE_CONT_H0_hypergeoR(12, 18, 2, 8, 4, 16)$qH0,
               expected = 0.005735633,
               tolerance = 1e-6)
})

test_that("Finding MLE for Compliers", {
  expect_equal(object = FindMLE_CONT_H1_hypergeoR(7, 4, 2, 3, 1, 5)$qH1,
               expected = 0.0857347,
               tolerance = 1e-6)
  expect_equal(object = FindMLE_CONT_H1_hypergeoR(16, 1, 5, 1, 2, 8)$qH1,
               expected = 0.06300415,
               tolerance = 1e-6)
  expect_equal(object = FindMLE_CONT_H1_hypergeoR(12, 18, 2, 8, 4, 16)$qH1,
               expected = 0.008749041,
               tolerance = 1e-6)
})