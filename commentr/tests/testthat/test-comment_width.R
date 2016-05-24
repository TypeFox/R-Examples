
context("comment_width")


test_that("Correct comment width", {
  expect_that(comment_width(), is_a("numeric"))
  expect_that(comment_width(45), is_equivalent_to(45))
  expect_that(comment_width(10000), is_less_than(200))
  expect_that(comment_width(10000), is_less_than(200))
  expect_that(comment_width("hej och hå"), throws_error())
  expect_that(comment_width("script_width"), equals(getOption("width") - 5))
  expect_that(comment_width("a4portrait"), is_a("numeric"))
  expect_that(comment_width("a4landscape"), is_a("numeric"))
  expect_that(comment_width("a3portrait"), is_a("numeric"))
  expect_that(comment_width("a3landscape"), is_a("numeric"))
  expect_that(comment_width("a4portrait", a4portrait_width = 40), equals(40))
})
