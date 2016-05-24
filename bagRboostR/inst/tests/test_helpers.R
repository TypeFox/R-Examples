context("Common helper functions")

test_that("outcomeLabel returns the formula outcome as a string", {
  formula <- as.formula("activity~.")
  expect_equal(outcomeLabel(formula),"activity")
  formula <- as.formula("activity ~ .")
  expect_equal(outcomeLabel(formula),"activity")
})