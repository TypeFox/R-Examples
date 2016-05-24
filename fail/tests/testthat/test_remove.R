context("remove and clear")

test_that("remove", {
  path = tempfile()
  f = fail(path)

  f$put(a = 1, b = 2)
  expect_equal(f$remove("b"), setNames(TRUE, "b"))
  expect_equal(f$ls(), "a")
  f$put(b = 2)
  expect_equal(f$ls(), letters[1:2])
  expect_equal(f$remove(letters[1:2]), setNames(c(TRUE, TRUE), letters[1:2]))
  f$put(a = 1, b = 2)

  # invalid keys and empty sets
  expect_warning(f$remove("c"), "Files not removed")
  expect_error(f$remove())
  expect_equal(f$remove(character(0L)), setNames(logical(0), character(0)))

  # cache
  expect_equal(f$get("a", use.cache=TRUE), 1)
  expect_equal(f$cached(), "a")
  expect_equal(f$remove("a"), setNames(TRUE, "a"))
  expect_equal(f$cached(), character(0L))
  f$put(a = 1)
  f$get("a", use.cache=TRUE)
  f$get("b", use.cache=TRUE)
  expect_equal(f$cached(), letters[1:2])
  f$remove(letters[1:2])
  expect_equal(f$cached(), character(0L))
})

test_that("clear", {
  path = tempfile()
  f = fail(path)
  f$put(a = 1, b = 2)
  f$get("a", use.cache=TRUE)
  f$get("b", use.cache=TRUE)
  expect_equal(f$cached(), letters[1:2])
  f$clear()
  expect_equal(f$cached(), character(0L))
})
