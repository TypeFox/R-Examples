library(testthat)

context("tests")

load("flo_results.RData")

test_that("metrics work as expected", {
  expect_less_than(sum(abs(influenceR::betweenness(flo_graph) - flo_bet)), 0.1) # some instability in betweenness results
  expect_equal(influenceR::ens(flo_graph), flo_ens)
  expect_equal(influenceR::constraint(flo_graph), flo_constraint)
  expect_equal(influenceR::bridging(flo_graph), flo_bridge)
})



