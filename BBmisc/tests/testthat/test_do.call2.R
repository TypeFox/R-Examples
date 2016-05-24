context("do.call2")

test_that("do.call2", {
  f = function(...) list(...)
  expect_equal(do.call2("f", a=1, b=2), list(a=1, b=2))
  expect_equal(do.call2("f", .args=list(a=1, b=2)), list(a=1, b=2))
  expect_equal(do.call2("f", a=1, .args=list(b=2)), list(a=1, b=2))

  df = iris
  expect_equal(do.call2("f", df), list(df))
  expect_equal(do.call2("f", .args = list(df)), list(df))

  f = function(x, data) data[[x]]
  expect_equal(do.call2("f", "Species", data=iris), iris$Species)
  expect_equal(do.call2("f", "Species", iris), iris$Species)
  expect_equal(do.call2("f", data = iris, "Species"), iris$Species)
  expect_equal(do.call2("f", "Species", .args = list(data = iris)), iris$Species)
  expect_equal(do.call2("f", data = iris, .args = list(x = "Species")), iris$Species)

  expect_error(do.call2(mean, 1:10), "string")
})
