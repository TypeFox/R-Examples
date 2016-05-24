context("Get Token")

token <- c(Sys.getenv('ClarifaiId'), Sys.getenv('ClarifaiSecret'))

test_that("get_token works well", {
  skip_on_cran()
  secret_id(token)
  get_token()
  expect_that(Sys.getenv("ClarifaiToken"), is_a("character"))
})