library(generator)
library(detector)
context("detect()")

test_df <-
  data.frame(phone_number = r_phone_numbers(1000),
             email_address = r_email_addresses(1000),
             national_identification_number = r_national_identification_numbers(1000),
             mixed = c(r_phone_numbers(300), r_email_addresses(200), r_national_identification_numbers(500)),
             stringsAsFactors = FALSE)

test_that("Produces the correct output.", {
})

test_that("Produces the correct output type.", {
  expect_is(detect(letters), "data.frame")
  expect_is(detect(as.Date("2014-01-01")), "data.frame")
  expect_is(detect(test_df), "data.frame")
})

test_that("Produces the correct messages, warnings, and errors.", {
})

