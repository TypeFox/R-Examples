context("mapValues")

test_that("mapValues", {
  expect_equal(mapValues(1:3, 2, 3), c(1, 3, 3))
  expect_equal(mapValues(letters[1:5], letters[1:5], rev(letters[1:5])), rev(letters[1:5]))
  expect_equal(mapValues(factor(c("a", "b", "c")), "b", "zzz"), factor(c("a", "zzz", "c"), levels = c("a", "zzz", "c")))
  expect_equal(mapValues(c("aab", "aba", "baa"), "aa", "zz", regex = TRUE), c("zzb", "aba", "bzz"))
  expect_equal(mapValues(c("aab", "aba", "baa"), "^aa.+", "zz", regex = TRUE), c("zz", "aba", "baa"))

  expect_error(mapValues(iris, 1, 1), "atomic")
  expect_error(mapValues(1:10, 1:2, 1), "length")
  expect_error(mapValues(1:10, 1, 1:2), "length")
})
