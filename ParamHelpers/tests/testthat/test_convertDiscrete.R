context("convertDiscrete")


test_that("discrete param NameToValue", {
  f = function(x) 2 * x
  p = makeDiscreteParam(id = "x", values = list(a = "char", b = 2L, c = 2.2, d = f, "e"))
  expect_equal(discreteNameToValue(p, "b"), 2L)
  expect_equal(discreteNameToValue(p, "c"), 2.2)
  expect_equal(discreteNameToValue(p, "e"), "e")
  expect_error(discreteNameToValue(p, ""), "Names not used")
  p = makeIntegerParam(id = "y")
  expect_error(discreteNameToValue(p, "a"))
})

test_that("discrete param ValueToName", {
  f = function(x) 2 * x
  p = makeDiscreteParam(id = "x", values = list(a = "char", b = 2L, c = 2.2, d = f, "e"))
  expect_equal(discreteValueToName(p, 2L), "b")
  expect_equal(discreteValueToName(p, 2.2), "c")
  expect_equal(discreteValueToName(p, "e"), "e")
  expect_equal(discreteValueToName(p, f), "d")
  expect_error(discreteValueToName(p, 3), "Value not found")
  expect_error(discreteNameToValue(p, 2))
  p = makeIntegerParam(id = "y")
  expect_error(discreteNameToValue(p, "a"))
})


test_that("discrete vec param NameToValue", {
  f = function(x) 2 * x
  p = makeDiscreteVectorParam(id = "x", len = 2, values = list(a = "char", b = 2L, c = c(2.2, 3.3), d = f, "e"))
  expect_equal(discreteNameToValue(p, c("a", "b")), list(a = "char", b = 2L))
  expect_equal(discreteNameToValue(p, c("e", "b")), list(e = "e", b = 2L))
  expect_equal(discreteNameToValue(p, c("c", "d")), list(c = c(2.2, 3.3), d = f))
})

test_that("discrete vec param ValueToName", {
  f = function(x) 2 * x
  p = makeDiscreteVectorParam(id = "x", len = 2, values = list(a = "char", b = 2L, c = 2.2, d = f, "e"))
  expect_equal(discreteValueToName(p, list(2L, "char")), c("b", "a"))
  expect_equal(discreteValueToName(p, list(2.2, f)), c("c", "d"))
})
