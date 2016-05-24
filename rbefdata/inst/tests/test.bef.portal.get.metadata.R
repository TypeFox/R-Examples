source("test.helpers.R")

context("Get metadata with bef.portal.get.metadata()")

test_that("it downloads metadata for john doe", {
  given_the_user_is("invalid")
  given_the_portal_is("development")
  expect_that(bef.portal.get.metadata(dataset_id=7), is_a("list"))
})

test_that("it throws an error if eml is not found", {
  given_the_portal_is("development")
  expect_error(bef.portal.get.metadata(dataset_id=1000), "failed to load HTTP resource")
})

test_that("it gets metadata by id", {
  given_the_portal_is("development")
  expect_that(bef.portal.get.metadata(dataset_id=7), is_a("list"))
})

test_that("it gets metadata by a full url", {
  given_the_portal_is("development")
  expect_that(bef.portal.get.metadata(dataset_id=7), is_a("list"))
})
