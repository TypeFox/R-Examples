context("Create S3 methods")

test_that("create_s3_method", {
  expect_is(create_s3_method(), "function")
  expect_is(create_s3_method(`-`), "function")
  expect_is(create_s3_method(`-.pin`), "function")
})