context("paramValueToString")

test_that("paramValueToString ", {
  u = makeNumericParam("u")
  v = makeIntegerVectorParam("v", len=2)
  w = makeDiscreteParam("w", values=list(a=1, b=list()))
  x = makeLogicalParam("x")
  y = makeLogicalVectorParam("y", len=2)
  z = makeDiscreteVectorParam("z", len=2, values=list(a=1, b=list()))
  s = makeCharacterParam("s")

  ps = makeParamSet(u, v, w, x, y, z, s)
  expect_equal(paramValueToString(u, 1), "1")
  expect_equal(paramValueToString(u, 1.2345), "1.23")
  expect_equal(paramValueToString(v, c(1, 2)), "1,2")
  expect_equal(paramValueToString(w, 1), "a")
  expect_equal(paramValueToString(w, list()), "b")
  expect_equal(paramValueToString(x, TRUE), "TRUE")
  expect_equal(paramValueToString(y, c(TRUE, FALSE)), "TRUE,FALSE")
  expect_equal(paramValueToString(z, list(1, 1)), "a,a")
  expect_equal(paramValueToString(z, list(1, list())), "a,b")
  expect_equal(paramValueToString(s, "PH"), "PH")
  expect_equal(paramValueToString(ps, list(u = 1, v = 1:2, w = list(), x = FALSE,
    y = c(TRUE, FALSE), z = list(1, list()), s = "PH")),
    "u=1; v=1,2; w=b; x=FALSE; y=TRUE,FALSE; z=a,b; s=PH")
})

test_that("requires works", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericParam("y", lower=1, upper=2, requires = quote(x == "a")),
    makeIntegerVectorParam("z", len=2, lower=1, upper=20, requires = quote(x == "b"))
  )
  expect_equal(paramValueToString(ps, list(x="a", y=1, z=NA), show.missing.values=TRUE), "x=a; y=1; z=NA")
  expect_equal(paramValueToString(ps, list(x="b", y=NA, z=2:3), show.missing.values=TRUE), "x=b; y=NA; z=2,3")
  expect_equal(paramValueToString(ps, list(x="b", z=2:3), show.missing.values=TRUE), "x=b; z=2,3")
  expect_equal(paramValueToString(ps, list(x="a", y=1, z=NA), show.missing.values=FALSE), "x=a; y=1")
  expect_equal(paramValueToString(ps, list(x="b", y=NA, z=2:3), show.missing.values=FALSE), "x=b; z=2,3")
  expect_equal(paramValueToString(ps, list(x="b", z=2:3), show.missing.values=FALSE), "x=b; z=2,3")
})

test_that("num.format works", {
  x = makeNumericParam("x")
  expect_equal(paramValueToString(x, 1), "1")
  expect_equal(paramValueToString(x, 1, num.format = "%.3f"), "1.000")
})

