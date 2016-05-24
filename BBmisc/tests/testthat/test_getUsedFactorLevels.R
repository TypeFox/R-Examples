context("getUsedFactorLevels")

test_that("getUsedFactorLevels", {
  used.levels = letters[1:3]
  all.levels = letters
  x = factor(x = used.levels, levels = all.levels)
  expect_equal(sort(getUsedFactorLevels(x)), used.levels)

  x = factor(x = used.levels)
  expect_equal(used.levels, levels(x))
})
