context("Test tag retrieval")

test_that("Basic tag queries can be executed", {
  skip_on_cran()
  results <- guardian_tags(query = "green", api_key = "test")
  expect_equal(length(results), 8)
  expect_equal(names(results), c("status", "userTier", "total", "startIndex", "pageSize", "currentPage", 
                                "pages", "results"))
})

test_that("More complex tag queries can be executed", {
  skip_on_cran()
  results <- guardian_tags(query = "green", section = "technology", api_key = "test")
  expect_equal(length(results), 8)
  expect_equal(names(results), c("status", "userTier", "total", "startIndex", "pageSize", "currentPage", 
                                 "pages", "results"))
})
