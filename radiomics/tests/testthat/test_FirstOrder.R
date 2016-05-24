library(radiomics)
context("First Order Statistics")

test_that("First order stats are correct", {
  expect_equal(calc_energy(hallbey), 42)
  expect_equal(calc_entropy(hallbey), 1.924, tolerance = .002)
  expect_equal(calc_kurtosis(hallbey), -1.416, tolerance = .002)
  expect_equal(calc_meanDeviation(hallbey), 0.90625, tolerance = .002)
  expect_equal(calc_skewness(hallbey), 0.1554, tolerance = .002)
  expect_equal(calc_uniformity(hallbey), 0.273, tolerance = .002)
})