context("comment")

test_that("comments", {
  expect_equal(format(make_comment("a")), "# a")
  expect_equal(format(make_comment("a", "b")), c("# a", "# b"))
  expect_equal(format(make_comment(c("a", "b"))), c("# a", "# b"))
})

test_that("error checking", {
  expect_error(make_comment(), "At least")
})
