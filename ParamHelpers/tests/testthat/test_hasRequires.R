context("hasRequires")

test_that("hasRequires", {
  ps = makeParamSet()
  expect_false(hasRequires(ps))
  ps = makeParamSet(
    makeNumericParam("x", lower=1, upper=2)
  )
  expect_false(hasRequires(ps))
  ps = makeParamSet(
    makeDiscreteParam("x", values=1:2),
    makeNumericParam("y", lower=1, upper=2, requires=quote(x == 1))
  )
  expect_true(hasRequires(ps))
})

