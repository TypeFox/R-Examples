context("apply")

test_that("apply", {
  path = tempfile()
  f = fail(path)

  f$put(a = 1:10, b = 1:100, c = 1:1000)
  x = f$apply(mean)
  expect_true(is.list(x))
  expect_true(all(names(x) %in% letters[1:3]))
  expect_equal(sort(unlist(x, use.names=FALSE)), c(5.5, 50.5, 500.5))

  # subsetting keys works
  x = f$apply(mean, keys=c("c", "b"))
  expect_equal(names(x), c("c", "b"))

  # simplify works
  x = f$apply(mean, keys=letters[1:3], simplify=TRUE)
  expect_equal(x, setNames(c(5.5, 50.5, 500.5), letters[1:3]))

  # use.names works
  x = f$apply(mean, keys=letters[1:3], use.names=FALSE)
  expect_true(is.null(names(x)))

  # passing arguments works
  fun = function(x, y) mean(x) + y
  x = f$apply(fun, y = -0.5, simplify=TRUE)
  expect_equal(x, setNames(c(5, 50, 500), letters[1:3]))

  # error handling
  f$remove(f$ls())
  f$put(a = 1, b = 2, c = "NA")
  expect_error(f$apply(log), "key 'c'")

  # invalid keys and empty sets
  expect_equal(length(f$apply(identity, keys=NULL)), 0)
  expect_error(f$apply(identity, keys="xxx"))
  expect_equal(length(f$apply(identity, keys=character(0L))), 0)
  f$remove(f$ls())
  expect_equal(length(f$apply(identity)), 0)
})
