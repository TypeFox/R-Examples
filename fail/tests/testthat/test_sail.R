context("sail")

test_that("sail", {
  # test basic usage here.
  # most of the stuff is identical with fail ...

  path = tempfile()
  s = sail(path)

  expect_equal(s$ls(), character(0L))
  expect_equal(s$put(a = 1, b = 2), letters[1:2])
  expect_equal(s$ls(), letters[1:2])
  expect_equal(s$put(li = list(c = 3)), letters[3])
  expect_equal(s$ls(), letters[1:3])
  s$put(d = 4, li = list(e = 5))
  expect_equal(s$ls(), letters[1:5])

  expect_equal(s$as.list(), setNames(as.list(1:5), letters[1:5]))

  s = sail(path, simplify = FALSE)
  x = s$as.list()
  expect_equal(x[[1]], list(a = 1))

  s = sail(path, simplify = TRUE)
  expect_equal(s$as.list(), s$apply(identity))

  expect_equal(s$cached(), character(0L))
  s$get("a", use.cache=TRUE)
  expect_equal(s$cached(), "a")
  expect_equal(s$get("a", use.cache = TRUE), 1)
  expect_equal(s$clear(), "a")
  expect_equal(s$cached(), character(0L))

  expect_equal(s$remove(s$ls()), setNames(rep(TRUE, 5L), letters[1:5]))
})
