source("test.helpers.R")

context("Get a dataset with bef.portal.get.dataset()")

test_that("it throws an error if the credentials are invalid", {
  given_the_user_is("invalid")
  given_the_portal_is("development")
  expect_error(bef.portal.get.dataset(dataset_id=7), "*not found or not accessible*")
})

test_that("it throws an error if the dataset is not found", {
  given_the_user_is("valid")
  given_the_portal_is("development")
  expect_error(bef.portal.get.dataset(dataset_id=1000), "*Internal Server Error*")
})

test_that("it gets a dataset by dataset_id", {
  given_the_user_is("valid")
  given_the_portal_is("development")
  expect_that(bef.portal.get.dataset(dataset_id=7), is_a("data.frame"))
})

test_that("it gets a dataset by full_url", {
  given_the_portal_is("development")
  expect_that(bef.portal.get.dataset(full_url= paste("http://befdatadevelepment.biow.uni-leipzig.de/datasets/7/download.csv?user_credentials=", given_the_user_is("valid"), sep = "")), is_a("data.frame"))
})

test_that("it gets multiple datasets by id", {
  given_the_user_is("valid")
  given_the_portal_is("development")
  expect_that(bef.portal.get.dataset(dataset_id=c(7,8)), is_a("list"))
})
