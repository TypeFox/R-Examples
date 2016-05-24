PrincipleBranchAnswers <- runif(2500, min = -1, max = 703.227)
PrincipleBranchTests <- PrincipleBranchAnswers * exp(PrincipleBranchAnswers)
SecondaryBranchAnswers <- runif(2500, min = -703.227, max = -1)
SecondaryBranchTests <- SecondaryBranchAnswers * exp(SecondaryBranchAnswers)

context("Testing lambertW")

test_that("Functions return proper values", {
  expect_equal(lambertW0(PrincipleBranchTests), PrincipleBranchAnswers)
  expect_equal(lambertWm1(SecondaryBranchTests), SecondaryBranchAnswers)
})

test_that("Function behaves properly near 0", {
  V0 <- seq(-1e-2, 1e-2, 1e-6)
  V0E <- V0 * exp(V0)
  LV0 <- lambertW0(V0E)
  expect_equal(V0, LV0)
  Vm1 <- V0[V0 < 0]
  LVm1 <- lambertWm1(Vm1)
  expect_equal(Vm1, LVm1 * exp(LVm1))
})

test_that("NaNs are returned for values outside domain", {
  expect_true(is.nan(lambertW0(-1)))
  expect_true(is.nan(lambertWm1(-1)))
  expect_equal(lambertWm1(0), -Inf)
  expect_true(is.nan(lambertWm1(1)))
  expect_true(is.nan(lambertW0(c(1, -1)))[[2]])
})