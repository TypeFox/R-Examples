context("Test queries")

test_that("basic project pageview queries work", {
  result <- project_pageviews()
  expect_true(is.data.frame(result))
  expect_true(nrow(result) == 1)
  expect_true(ncol(result) == 6)
})

test_that("Basic top-article queries work", {
  result <- top_articles()
  expect_true(is.data.frame(result))
  expect_true(nrow(result) == 1000)
  expect_true(ncol(result) == 8)
})

test_that("Basic per-article queries work", {
  result <- article_pageviews()
  expect_true(is.data.frame(result))
  expect_true(nrow(result) == 1)
  expect_true(ncol(result) == 6)
})
