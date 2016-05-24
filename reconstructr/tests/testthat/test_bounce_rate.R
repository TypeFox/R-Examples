context("bounce_rate")

test_that("bounce_rate handles 0-entry lists elegantly", {
  expect_that(is.nan(bounce_rate(list())), equals(TRUE))
})

test_that("bounce_rate handles 1-entry, all-bounce lists elegantly", {
  expect_that(bounce_rate(list(12)), equals(100))
})

test_that("bounce_rate handles 1-entry, no-bounce lists elegantly", {
  expect_that(bounce_rate(list(c(12,300))), equals(0))
})

test_that("bounce_rate handles multi-entry, all-bounce lists elegantly", {
  expect_that(bounce_rate(list(12,300)), equals(100))
})

test_that("bounce_rate handles multi-entry, no-bounce lists elegantly", {
  expect_that(bounce_rate(list(c(12,300),c(30,90))), equals(0))
})

test_that("bounce_rate handles multi-entry, mixed-bounce lists elegantly", {
  expect_that(bounce_rate(list(c(12,300),30)), equals(50))
})