context("lsort")

test_that("lsort", {
  expect_equal(lsort(c("c", "a", "b")), c("a", "b", "c"))
  expect_equal(lsort( c("a", "ä", "ö", "o")),  c("a", "o", "ä", "ö"))
})
