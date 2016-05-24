library(junr)
library(httr)

context("Get data")

base_url <- "http://api.datosabiertos.presidencia.go.cr/api/v2/datastreams/"
api_key <- "0bd55e858409eefabc629b28b2e7916361ef20ff"

test_that("The connection to the test url gets a response", {
  r <- GET(paste(base_url, "?auth_key=", api_key, sep = ""), accept_json())
  expect_true(r$status_code %in% c(200, 403, 500))
})

test_that("The data index is read correctly", {
  test_index <- get_index(base_url, api_key)
  expect_true(exists("test_index"))
})

