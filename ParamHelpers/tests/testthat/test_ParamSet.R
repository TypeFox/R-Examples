context("ParameterSet")

test_that("empty paramset", {
  ps = makeParamSet()
  expect_equal(getParamIds(ps), character(0))
  expect_equal(getParamLengths(ps), integer(0))
  expect_equal(getLower(ps), numeric(0))
  expect_equal(getUpper(ps), numeric(0))
  expect_equal(getValues(ps), list())
  expect_output(print(ps), "Empty parameter set.")
  expect_true(isEmpty(ps))
})


test_that("mixed paramset", {
  ps = makeParamSet(
    makeNumericParam("u", lower = 1),
    makeIntegerParam("v", lower = 1, upper = 2),
    makeDiscreteParam("w", values = 1:2),
    makeLogicalParam("x"),
    makeDiscreteVectorParam("y", len = 2, values = c("a", "b"))
  )
  expect_equal(getParamIds(ps), c("u", "v", "w", "x", "y"))
  expect_equal(getParamLengths(ps), c(u = 1L, v = 1L, w = 1L, x = 1L, y = 2L))
  expect_equal(getLower(ps), c(u = 1, v = 1L))
  expect_equal(getUpper(ps), c(u = Inf, v = 2L))
  values1 = list(1,2); names(values1) = 1:2
  values2 = list("a","b"); names(values2) = c("a", "b")
  values3 = list(TRUE, FALSE); names(values3) = c("TRUE", "FALSE")
  expect_equal(getValues(ps), list(w = values1, x = values3, y = values2))
  expect_output(print(ps), "u\\s*numeric")
})

test_that("mixed paramset 2", {
  p1 = makeNumericParam(id = "x1", lower = -1L, upper = 1)
  p2 = makeNumericParam(id = "x2", lower = 0, upper = Inf)
  p3 = makeDiscreteParam(id = "x3", values = list(a = "char", b = 2L, c = 2.2, "e"))
  ps = makeParamSet(p1, p2, p3)

  expect_true(isFeasible(ps, list(0,0,"char")))
  expect_true(!isFeasible(ps, list(2,0,"char")))
  expect_equal(getLower(ps), c(x1 = -1, x2 = 0))
  expect_equal(getUpper(ps), c(x1= 1, x2 = Inf))
  expect_true(isFeasible(ps, list(x3 = 2L, x1 = 0)))
})

test_that("cannot build param set from wrong stuff", {
  expect_error(makeParamSet(1), "Param")
  expect_error(
    makeParamSet(
      makeNumericParam("x"),
      makeNumericParam("x")
    ), "unique"
 )
})

test_that("paramset with vector", {
  ps = makeParamSet(
    makeNumericParam("x", lower = 1),
    makeIntegerVectorParam("y", len = 2, lower = 3L, upper = 9L)
  )
  expect_equal(getParamIds(ps), c("x", "y"))
  expect_equal(getParamLengths(ps), c(x = 1L, y = 2L))
  expect_equal(getParamIds(ps, repeated = TRUE, with.nr = FALSE), c("x", "y", "y"))
  expect_equal(getParamIds(ps, repeated = TRUE, with.nr = TRUE), c("x", "y1", "y2"))
  expect_equal(getLower(ps), c(x = 1, y = 3L, y = 3L))
  expect_equal(getUpper(ps), c(x = Inf, y = 9L, y = 9L))
  expect_true(is.list(getValues(ps)))
  expect_equal(length(getValues(ps)), 0)
  expect_output(print(ps), "x\\s*numeric")
})

test_that("list works", {
  ps = makeParamSet(params = list(
    makeNumericParam("u", lower = 1),
    makeIntegerParam("v", lower = 1, upper = 2)
  ))
  expect_equal(getParamIds(ps), c("u", "v"))
  expect_output(print(ps), "u\\s*numeric")
})

test_that("combination with c works", {
  ps1 = makeParamSet(
    makeNumericParam("u", lower = 1),
    makeIntegerVectorParam("v", len = 2, lower = 3L, upper = 9L)
  )
  ps2 = makeParamSet(
    makeDiscreteParam("w", values = 1:2),
    makeLogicalParam("x")
  )
  # ps with "difficult" name
  ps3 = makeParamSet(
    makeNumericParam("params")
  )
  ps = c(ps1, ps2)
  expect_equal(length(ps$pars), 4)
  expect_equal(getParamIds(ps), c(getParamIds(ps1), getParamIds(ps2)))

  ps = c(ps1, ps3)
  expect_equal(length(ps$pars), 3)
  expect_equal(getParamIds(ps), c(getParamIds(ps1), getParamIds(ps3)))
})


test_that("cannot combine with overlapping names", {
  ps1 = makeParamSet(
    makeNumericParam("u")
  )
  ps2 = makeParamSet(
    makeDiscreteParam("u", values = 1:2)
  )
  expect_error(c(ps1, ps2), "unique")
})


test_that("unknown length works", {
  ps = makeParamSet(
    makeNumericLearnerParam("u", lower = 2),
    makeNumericVectorLearnerParam("v", len = NA, lower = 1)
  )
  expect_equal(getParamLengths(ps), c(u = 1, v = NA))
  expect_true(isFeasible(ps, list(3, c(2))))
  expect_true(isFeasible(ps, list(3, c(2, 4))))
  expect_false(isFeasible(ps, list(3, c(2, 0))))
  expect_error(sampleValue(ps), "Cannot sample")
})


test_that("makeNumericParamset", {
  ps = makeNumericParamSet(lower = 10, len = 2, vector = TRUE)
  expect_equal(getParamIds(ps), "x")
  expect_equal(getParamLengths(ps), c(x = 2))
  expect_equal(getLower(ps), c(x = 10, x = 10))
  expect_equal(getUpper(ps), c(x = Inf, x = Inf))
  ps = makeNumericParamSet(id = "y", upper = 3, len = 2, vector = FALSE)
  expect_equal(getParamIds(ps), c("y1", "y2"))
  expect_equal(getParamLengths(ps), c(y1 = 1, y2 = 1))
  expect_equal(getLower(ps), c(y1 = -Inf, y2 = -Inf))
  expect_equal(getUpper(ps), c(y1 = 3, y2 = 3))
  ps = makeNumericParamSet(lower = c(1,2), vector = TRUE)
  expect_equal(getLower(ps), c(x = 1, x = 2))
  ps = makeNumericParamSet(lower = c(1,4), upper = c(2,7), vector = FALSE)
  expect_equal(getLower(ps), c(x1 = 1, x2 = 4))
  expect_equal(getUpper(ps), c(x1 = 2, x2 = 7))
})

test_that("requires works", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b"), default = "a"),
    makeNumericParam("y", requires = quote(x == "a")),
    makeIntegerVectorParam("z", len = 2, requires = quote(x == "b"))
  )
  expect_true(isFeasible(ps, list(x = "a", y = 1,  z = NA)))
  expect_false(isFeasible(ps, list(x = "a", y = NA, z = 1)))
  expect_false(isFeasible(ps, list(x = "a", y = NA, z = c(NA, NA))))
  expect_false(isFeasible(ps, list(x = "b", y = 1, z = c(2,2))))
  expect_true(isFeasible(ps, list(x = "b", y = NA, z = c(2,2))))
  expect_true(isRequiresOk(ps, list(x = "a")))
  expect_true(isRequiresOk(ps, list(x = "c"))) #out of bound is supposed to be ignored
  expect_true(isRequiresOk(ps, list(y = 1)))
  expect_error(isRequiresOk(ps, list(y = 1), use.defaults = FALSE))
  expect_error(isRequiresOk(ps, list(y = 1, x = "b")), 'x == "a"')
})

test_that("requires chains work", {
  ps = makeParamSet(
    makeLogicalLearnerParam("a", default = FALSE),
    makeLogicalLearnerParam("b", default = FALSE, requires = quote(a == TRUE)),
    makeLogicalLearnerParam("c", default = FALSE, requires = quote(b == TRUE))
  )
  expect_true(isFeasible(ps, list(a = FALSE, b = NA, c = NA)))
  expect_true(isFeasible(ps, list(a = TRUE, b = FALSE, c = NA)))
  expect_true(isFeasible(ps, list(a = TRUE, b = TRUE, c = FALSE)))
  expect_true(isFeasible(ps, list(a = TRUE, b = TRUE, c = TRUE)))
})

test_that("print works", {
  ps = makeParamSet(
    makeIntegerParam("ntree", lower = 10, upper = 50),
    makeNumericVectorParam("cutoff", len = 3, lower = 0.001, upper = 1, trafo = function(x) 0.9*x/sum(x))
  )
  expect_output(print(ps), "numericvector")
  expect_output(print(ps), "3")

  ps = makeParamSet(
    makeIntegerLearnerParam(id = "x", default = 50L, lower = 1L),
    makeNumericVectorLearnerParam(id = "v", default = 1:2, len = 2L)
  )
  expect_output(print(ps, trafo = FALSE, used = FALSE), "numericvector")
})

