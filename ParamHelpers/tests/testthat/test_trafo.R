context("trafo")

test_that("trafoValue with param", {
  p = makeNumericParam(id="x", lower=-10, upper=10, trafo=function(x) x^2)
  expect_equal(trafoValue(p, 1), 1)
  expect_equal(trafoValue(p, -5), 25)
})

test_that("trafoValue with param set", {
  ps = makeParamSet(
    makeIntegerParam("u", trafo=function(x) 2*x),
    makeNumericVectorParam("v", len=2, trafo=function(x) x/sum(x)),
    makeDiscreteParam("w", values=c("a", "b"))
  )
  expect_equal(trafoValue(ps, list(3, c(2, 4), "a")), list(u=6, v=c(2/6, 4/6), w="a"))
})


test_that("trafo opt.path", {
  ps = makeParamSet(
    makeNumericParam("x", lower=-2, upper=2, trafo=function(x) 2^x)
  )
  op = makeOptPathDF(par.set = ps, y.names = "y", minimize = TRUE, add.transformed.x = FALSE)
  addOptPathEl(op, x=list(x = -2), y = 0)
  addOptPathEl(op, x=list(x = 2), y = 0)
  expect_error(addOptPathEl(op, x = list(x = 3), y = 0), "infeasible")
  op2 = trafoOptPath(op)
  df = as.data.frame(op2)
  expect_equal(df$x, c(1/4, 4))
  expect_error(trafoOptPath(op2), "Cannot further trafo")

  ps = makeParamSet(
    makeIntegerParam("u", trafo=function(x) 2*x),
    makeNumericVectorParam("v", len = 2, trafo=function(x) x/sum(x)),
    makeDiscreteParam("w", values=c("a", "b"))
  )
  op = makeOptPathDF(ps, "y", TRUE)
  addOptPathEl(op, x=list(3, c(2, 4), "a"), y = 0, dob = 1, eol = 1)
  addOptPathEl(op, x=list(4, c(5, 3), "b"), y = 2, dob = 5, eol = 7)
  op2 = trafoOptPath(op)
  df2 = as.data.frame(op2)
  df2b = rbind(
    data.frame(u = 6, v1 = 2/6, v2 = 4/6, w = "a", y = 0, dob = 1, eol = 1,
      stringsAsFactors = TRUE),
    data.frame(u = 8, v1 = 5/8, v2 = 3/8, w = "b", y = 2, dob = 5, eol = 7,
      stringsAsFactors = TRUE)
  )
  expect_equal(df2, df2b)
})

