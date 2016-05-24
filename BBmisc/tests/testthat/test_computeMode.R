context("computeMode")

test_that("computeMode", {
  # factor
  expect_equal(computeMode(as.factor(c(1:2, 2L, 2L))), "2")
  #character
  expect_equal(computeMode(c("1","2","3"), ties.method="last"), "3")
  # numeric
  expect_equal(computeMode(c(1,1,2,3)), 1)
  # integer
  expect_equal(computeMode(c(1:2, 2L, 2L), ties.method="first"), 2L)
  expect_equal(computeMode(c(1:2, 2L, 2L), ties.method="random"), 2L)
  expect_equal(computeMode(c(1:2, 2L, 2L), ties.method="last"), 2L)
  # logical
  expect_equal(computeMode(c(TRUE, FALSE, FALSE)), FALSE)
  expect_equal(computeMode(c(TRUE, TRUE, FALSE)), TRUE)

  # na.rm
  expect_equal(computeMode(c(1,1,2,3, NA, NA, NA), na.rm=FALSE), as.numeric(NA))
  expect_equal(computeMode(c(1,1,2,3, NA, NA, NA), na.rm=TRUE), 1)
})
