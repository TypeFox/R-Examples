
context("Utilities")

test_that("%||%", {

  expect_equal(1 %||% "what", 1)
  expect_equal(NULL %||% "what", "what")
  expect_equal(numeric() %||% "what", numeric())

})

test_that("null_or_any_na", {

  expect_true(null_or_any_na(NULL))
  expect_true(null_or_any_na(c(1, NA, 2)))

  expect_false(null_or_any_na(list()))
  expect_false(null_or_any_na(numeric()))
  expect_false(null_or_any_na("foo"))
})

test_that("drop_last", {

  expect_equal(drop_last(1:10), 1:9)
  expect_equal(drop_last(1), numeric())
  expect_equal(drop_last(list(1)), list())
  expect_equal(drop_last(list(1, 2)), list(1))
})
