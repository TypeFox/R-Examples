context("getParamNr")

test_that("getParamNr", {
  ps = makeParamSet()
  expect_equal(getParamNr(ps), 0L)
  expect_equal(getParamNr(ps, devectorize = TRUE), 0L)

  ps = makeParamSet(
    makeNumericParam("u", lower = 1),
    makeIntegerParam("v", lower = 1, upper = 2),
    makeDiscreteParam("w", values = 1:2),
    makeLogicalParam("x"),
    makeDiscreteVectorParam("y", len = 2, values = c("a", "b"))
  )
  expect_equal(getParamNr(ps), 5L)
  expect_equal(getParamNr(ps, devectorize = TRUE), 6L)
})


