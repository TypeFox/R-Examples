context("data can be loaded")
test_that("data can be loaded", {
  data(chelonia.simmap)
  expect_that(length(chelonia.simmap$tree$tip.label),equals(226))
})
