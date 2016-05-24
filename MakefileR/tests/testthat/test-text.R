context("text")

test_that("text", {
  expect_equal(format(make_text("a")), "a")
  expect_equal(format(make_text("a", "b")), c("a", "b"))
  expect_equal(format(make_text(letters)), letters)
})

test_that("error checking", {
  expect_error(make_text(), "At least")
})
