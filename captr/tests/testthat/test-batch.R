context("Create Batch")

token <- Sys.getenv('CaptricityToken')

test_that("batch creation happens successfully", {
  skip_on_cran()
  set_token(token)
  batch <- create_batch("new_batch")
  expect_that(batch, is_a("list"))
})