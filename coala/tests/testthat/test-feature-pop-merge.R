context("Feature Pop Merge")

test_that("Creation of merge features works", {
  feat <- feat_pop_merge(2, 2, 1)
  expect_equal(feat$get_time(), "par(2)")
  expect_equal(feat$get_population(), c(from = 2, to = 1))
})


test_that("generating scrm cmd for pop merge works", {
  model <- coal_model(4:5, 1) + feat_pop_merge(par_range("tau", 1, 2), 2, 1)
  expect_equal(get_simulator("scrm")$get_cmd(model),
               "scrm 9 1 -I 2 4 5 -ej tau 2 1 ")
})
