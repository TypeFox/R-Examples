source("test.helpers.R")

context("Get files associated to a dataset with bef.portal.get.attachment()")

test_that("it throws an error if there is no associated files", {
  given_the_portal_is(environment = "development")
  given_the_user_is(condition = "valid")
  id = given_the_dataset_is(available = FALSE)
  expect_error(bef.portal.get.attachments(dataset = id), "*Internal Server Error*")
})

test_that("it creates a new folder for the download files", {
  given_the_portal_is(environment = "development")
  given_the_user_is(condition = "valid")
  folder = bef.options("download_dir")
  id = given_the_dataset_is(available = TRUE)
  bef.portal.get.attachments(dataset = id)
  expect_that(file.exists(folder), is_true())
  if (file.exists(folder)) {
    unlink(folder, recursive = TRUE)
  }
})
