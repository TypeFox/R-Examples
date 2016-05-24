context("download successful")

test_that("basic download", {
  expect_is(webuse("auto"), "data.frame")
})

test_that("download version", {
  expect_is(webuse("auto", version = 14), "data.frame")
  expect_is(webuse("auto", version = 8), "data.frame")
})

test_that("data.frame added to .GlobalEnv", {
  webuse("lifeexp")
  expect_true("lifeexp" %in% ls())
})

test_that("data.frame added to specified environment", {
  e <- new.env()
  webuse("lifeexp", envir = e)
  expect_true("lifeexp" %in% ls(e))
})
