
context("Utilities")

test_that("merge_named_lists", {

  cases <- list(
    list(
      list(a = "a", b = "b"),
      list(a = "x", b = "y"),
      list(a = c("a", "x"), b = c("b", "y"))
    ),
    list(
      list(a = "a", b = "b"),
      list(),
      list(a = "a", b = "b")
    ),
    list(
      list(),
      list(a = "a", b = "b"),
      list(a = "a", b = "b")
    ),
    list(
      list(a = "a", b = "b"),
      list(c = "c"),
      list(a = "a", b = "b", c = "c")
    )
  )

  for (e in cases) {
    expect_equal(
      merge_named_lists(e[[1]], e[[2]]),
      e[[3]]
    )
  }
})

test_that("%||%", {

  expect_equal(
    "foo" %||% "bar",
    "foo"
  )
  expect_equal(
    NULL %||% "bar",
    "bar"
  )
  expect_equal(
    NULL %||% NULL,
    NULL
  )
})
