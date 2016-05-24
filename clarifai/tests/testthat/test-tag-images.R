context("Tag Local Image")

token <- c(Sys.getenv('ClarifaiId'), Sys.getenv('ClarifaiSecret'))

image <- system.file("extdata/metro-north.jpg", package = "clarifai")

test_that("tag_images works ok", {
  skip_on_cran()
  secret_id(token)
  get_token()
  tag <- tag_images(image)
  expect_that(tag, is_a("data.frame"))
})