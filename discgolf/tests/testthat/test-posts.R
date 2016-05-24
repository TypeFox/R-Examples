context("posts")

test_that("fails well with no input", {
  skip_on_cran()

  aa <- post_get(90)
  bb <- post_get(120)
  cc <- post_get(130)

  expect_is(aa, "list")
  expect_is(aa$name, "character")
  expect_is(aa$actions_summary, "data.frame")
  expect_equal(aa$username, "sckott")

  expect_is(bb, "list")
  expect_is(bb$name, "character")
  expect_is(bb$actions_summary, "data.frame")
  expect_equal(bb$username, "bw4sz")

  expect_is(cc, "list")
  expect_is(cc$name, "character")
  expect_is(cc$actions_summary, "data.frame")
  expect_equal(cc$username, "sckott")
})

test_that("fails well with no input", {
  skip_on_cran()

  expect_error(post_get(), "argument \"post_id\" is missing")
})

test_that("fails well with non-existent page", {
  skip_on_cran()

  expect_error(post_get("asfafsfadfasdfd"),
               "404 - The requested URL or resource could not be found.")
})

test_that("httr curl options work", {
  skip_on_cran()

  library("httr")
  expect_error(post_get("asdfadf", config = timeout(seconds = 0.001)))
})
