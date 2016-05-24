context("list, get, put")

test_that("list, get, put", {
  path = tempfile()
  f = fail(path)

  expect_equal(f$ls(), character(0L))
  expect_equal(f$put(a = 1, b = 2), letters[1:2])
  expect_equal(f$ls(), letters[1:2])
  expect_equal(f$put(li = list(c = 3)), letters[3])
  expect_equal(f$ls(), letters[1:3])
  f$put(d = 4, li = list(e = 5))
  expect_equal(f$ls(), letters[1:5])

  expect_equal(f$get("a"), 1)
  x = f$as.list()
  y = setNames(as.list(1:5), letters[1:5])
  expect_equal(x, y)

  # positional arguments
  expect_equal(f$pos(), 1)
  expect_equal(f$pos(2), 2)
  expect_equal(f$pos(6), NULL)

  path = tempfile()
  f = fail(path)
  f$put(1, 2, 3, keys = c("x", "y", "z"))
  expect_equal(f$ls(), c("x", "y", "z"))
  f$remove(f$ls())
  f$put(1, 2, 3, li = list(foo = 5), keys = c("x", "y", "z"))
  expect_equal(f$get("x"), 1)
  expect_equal(f$get("foo"), 5)

  # pattern works
  expect_equal(f$ls("^[xy]"), c("x", "y"))
  expect_equal(f$ls("a"), character(0L))

  # invalid keys and empty sets
  expect_error(f$get())
  expect_error(f$get("not_existing"))
  expect_error(f$put(li=list("a - b" = 1)))
  expect_equal(f$put(), character(0L))

  # cache
  expect_equal(f$cached(), character(0L))
  f$get("x", use.cache=TRUE)
  expect_equal(f$cached(), "x")
  f$clear()
  expect_equal(f$cached(), character(0L))
})
