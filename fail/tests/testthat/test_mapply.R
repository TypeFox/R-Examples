context("mapply")

test_that("mapply", {
  path = tempfile()
  f = fail(path)

  f$put(a = 1:10, b = 1:100)
  x = f$mapply(function(key, value) mean(value))
  expect_equal(x, list(a = 5.5, b = 50.5))
  x = f$mapply(function(key, value) mean(value), keys = c("b", "a"))
  expect_equal(x, list(b = 50.5, a = 5.5))
  x = f$mapply(function(key, value) mean(value), use.names=FALSE)
  expect_equal(x, list(5.5, 50.5))
  x = f$mapply(function(key, value) mean(value), simplify=TRUE)
  expect_equal(x, setNames(c(5.5, 50.5), letters[1:2]))
  x = f$mapply(function(key, value, y) mean(value - y), y = 1)
  expect_equal(x, list(a = 4.5, b = 49.5))
  x = f$mapply(function(key, value, y) mean(value - y), moreArgs = list(y = 1))
  expect_equal(x, list(a = 4.5, b = 49.5))

  expect_error(f$mapply(function(key, value, y) mean(value - y), y = 0, moreArgs = list(y = 1)))

  # error handling
  f$remove(f$ls())
  f$put(a = 1, b = 2, c = "NA")
  expect_error(f$mapply(function(key, value) log(value)), "key 'c'")


  # invalid keys and empty sets
  fun = function(key, value) value
  expect_equal(length(f$mapply(fun, keys=NULL)), 0)
  expect_error(f$mapply(fun, keys="xxx"))
  expect_equal(length(f$mapply(fun, keys=character(0L))), 0)
  f$remove(f$ls())
  expect_equal(length(f$mapply(fun)), 0)
})
