context("dropParams")

test_that("dropParams", {
  ps = makeParamSet(
    makeNumericParam("u", lower = 1),
    makeIntegerParam("v", lower = 1, upper = 2),
    makeDiscreteParam("w", values = 1:2)
  )
  expect_equal(getParamIds(dropParams(ps, "u")), c("v", "w"))
  expect_equal(getParamIds(dropParams(ps, c("u","v"))), "w")
  expect_output(dropParams(ps, c("u", "v", "w")), "Empty parameter set.")
})
