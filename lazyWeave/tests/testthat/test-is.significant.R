context("is.significant")

test_that("test significance at 0.05",
{
  expect_equal(
    is_significant(c(.9, .8, .65, .33, .10, .06, .051, .05, .049, .001), .05),
    c(rep(FALSE, 7), rep(TRUE, 3)))
})

test_that("test significance at 0.05",
{
  expect_equal(
    is_significant(c(.9, .8, .65, .33, .10, .06, .051, .05, .049, .001), .10),
    c(rep(FALSE, 4), rep(TRUE, 6)))
})