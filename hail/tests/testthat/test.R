context("hail")

test_that("hail_hydra works", {
  data <- hail_hydra("Walmart Eco Roof")
  expect_equal(ncol(data), 27)
  expect_equal(class(data), "data.frame")

})
