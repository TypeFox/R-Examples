context("coalesce")

test_that("coalesce", {
  expect_identical(coalesce(NULL), NULL)
  expect_identical(coalesce(1, NULL), 1)
  expect_identical(coalesce(NULL, 1), 1)
  f = function(a, b, c) coalesce(a,b,c)
  expect_identical(f(b = NULL, c = 1), 1)
})

test_that("coalesce works with functions", {
  ff = function(x) 123
  expect_identical(coalesce(NULL, ff), ff)
  # FIXME: this is a problem in R / missing. see issue 48
  #expect_identical(coalesce(NULL, min), min)
})


