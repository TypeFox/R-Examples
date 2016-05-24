context("getDefaults")

test_that("getDefaults", {
  ps = makeParamSet()
  expect_equal(getDefaults(ps), list())

  ps = makeParamSet(
    makeNumericParam("u", default = 2),
    makeIntegerParam("v"),
    makeCharacterParam("s", default = "foo"),
    makeUntypedParam("w", default = iris)
  )
  expect_equal(getDefaults(ps), list(u = 2, s = "foo", w = iris))

  ps = makeParamSet(
    makeIntegerParam("v")
  )
  expect_equal(getDefaults(ps), list())

  # test include.null
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b"), default = "a"),
    makeNumericVectorParam("y", len = 2),
    makeIntegerParam("z", default = 99)
  )
  expect_equal(
    getDefaults(ps, include.null = TRUE),
    list(x = "a", y = NULL, z = 99)
  )
  expect_equal(
    getDefaults(ps, include.null = FALSE),
    list(x = "a", z = 99)
  )
})

