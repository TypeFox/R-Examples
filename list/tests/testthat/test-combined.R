context("Tests combined estimator")
rm(list=ls())

test_that("combinedListDirect works", {
  data("combinedListExps")
  combinedListExps <- na.omit(combinedListExps)
  out.1 <- combinedListDirect(list1N ~ list1treat, 
                              data = subset(combinedListExps, directsfirst==1), 
                              treat = "list1treat", direct = "direct1")
  expect_that(round(out.1$comb.est, 8), is_equivalent_to(0.66561842))
  expect_that(round(out.1$placebo.I$estimate, 6), is_equivalent_to(1.053962))
  expect_that(round(out.1$placebo.II$estimate, 8), is_equivalent_to(0.06626651))
})
