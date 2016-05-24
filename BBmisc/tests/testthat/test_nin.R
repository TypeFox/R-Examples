context("nin")

test_that("nin", {
  expect_true(1 %nin% 2:3)
  expect_false(1 %nin% 1)
  expect_false(1 %nin% c(NA, 1))
  expect_true(1 %nin% c(NA, 2))
  expect_false(NA %nin% c(NA, 1))
  expect_true(NA %nin% 1:2)
})
