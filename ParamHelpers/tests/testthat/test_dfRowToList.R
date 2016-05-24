context("dfRowToList")


test_that("dfRowToList", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = list(a = iris, b = 123)),
    makeNumericParam("y", lower = 1, upper = 2),
    makeNumericVectorParam("z", len = 2, lower = 10, upper = 20)
  )
  des = generateDesign(10, par.set = ps)
  vals = dfRowsToList(des, ps)
  expect_true(is.list(vals) && length(vals) == nrow(des))
  for (i in seq_row(des)) {
    v = vals[[i]]
    expect_true(is.list(v) && length(v) == 3 && names(v) == c("x", "y", "z"))
    expect_true(identical(v$x, iris) || identical(v$x, 123))
    expect_true(is.numeric(v$y) && length(v$y) == 1 && v$y >= 1 && v$y <= 2)
    expect_true(is.numeric(v$z) && length(v$z) == 2 && all(v$z >= 10 & v$z <= 20))
  }
})

test_that("requires works", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericParam("y", lower = 1, upper=2, requires = quote(x == "a")),
    makeNumericVectorParam("z", len = 2, lower = 10, upper = 20, requires = quote(x == "b"))
  )
  des = generateDesign(10, par.set=ps)

  mycheck = function(vals) {
    expect_true(is.list(vals) && length(vals) == nrow(des))
    for (i in seq_row(des)) {
      v = vals[[i]]
      expect_true(is.list(v) && length(v) == 3L && (names(v) %in% c("x", "y", "z")))
      expect_true(is.character(v$x) && length(v$x) == 1 && v$x %in% c("a", "b"))
      if (v$x == "a") {
        expect_true(is.numeric(v$y) && length(v$y) == 1 && v$y >= 1 && v$y <= 2)
        expect_true(isScalarNA(v$z))
      } else if (v$x == "b") {
        expect_true(isScalarNA(v$y))
        expect_true(is.numeric(v$z) && length(v$z) == 2 && all(v$z >= 10 & v$z <= 20))
      }
    }
  }
  vals = dfRowsToList(des, ps)
  mycheck(vals)
})

test_that("ints in data frame work", {
  df = data.frame(x = 1:3, y = as.numeric(1:3))
  ps = makeParamSet(
    makeIntegerParam("x"),
    makeIntegerParam("y")
  )
  xs = dfRowsToList(df, ps)
  expect_equal(xs[[1]], list(x = 1L, y = 1L))
})

test_that("enforce.col.types works", {
  df = data.frame(u = c(NA, NA), x = 1:2, y = as.numeric(1:2), z = c("TRUE", "FALSE"))
  ps = makeParamSet(
    makeNumericParam("u"),
    makeIntegerParam("x"),
    makeIntegerParam("y"),
    makeLogicalParam("z")
  )
  x = dfRowToList(df, ps, 1L, enforce.col.types = TRUE)
  expect_equal(x, list(u = NA, x = 1L, y = 1L, z = TRUE))
})



