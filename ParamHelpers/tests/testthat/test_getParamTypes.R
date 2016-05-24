context("getParamTypes")

test_that("getParamTypes", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericVectorParam("y", len = 2),
    makeIntegerParam("z"),
    makeCharacterVectorParam("s", len = 2)
  )

  expect_equal(
    getParamTypes(ps, df.cols = FALSE, use.names = FALSE),
    c("discrete", "numericvector", "integer", "charactervector")
  )
  expect_equal(
    getParamTypes(ps, df.cols = TRUE, use.names = FALSE),
    c("factor", "numeric", "numeric", "integer", "character", "character")
  )
  expect_equal(
    getParamTypes(ps, df.cols = FALSE, use.names = TRUE),
    c(x = "discrete", y = "numericvector", z = "integer", s = "charactervector")
  )
  expect_equal(
    getParamTypes(ps, df.cols = TRUE, use.names = TRUE, with.nr = FALSE),
    c(x = "factor", y = "numeric", y = "numeric", z = "integer", s = "character",
      s = "character")
  )
  expect_equal(
    getParamTypes(ps, df.cols = TRUE, use.names = TRUE, with.nr = TRUE),
    c(x = "factor", y1 = "numeric", y2 = "numeric", z = "integer",
      s1 = "character", s2 = "character")
  )
})
