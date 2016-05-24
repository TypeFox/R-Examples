context("Get App. Info.")

token <- c(Sys.getenv('ClarifaiId'), Sys.getenv('ClarifaiSecret'))

test_that("get_info happens successfully", {
  skip_on_cran()
  secret_id(token)
  get_token()
  get_info <- get_info()
  expect_that(get_info, is_a("list"))
})