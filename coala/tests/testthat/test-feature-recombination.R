context("Feature Recombination")

test_that("Creation of growth features works", {
  expect_equal(feat_recombination(2)$get_rate(), "par(2)")
})


test_that("generating scrm cmd for growth works", {
  scrm <- get_simulator("scrm")
  model <- coal_model(15, 1, 25) + feat_recombination(par_range("rho", 0, 1))
  expect_equal(scrm$get_cmd(model), "scrm 15 1 -r rho locus_length ")
  model <- coal_model(15, 1, 25) + feat_recombination(5)
  expect_equal(scrm$get_cmd(model), "scrm 15 1 -r 5 locus_length ")
})


test_that("simulating recombination works", {
  model <- coal_model(5, 1, 10) + feat_recombination(1) + feat_mutation(1)
  expect_that(simulate(model), is_a("list"))
})
