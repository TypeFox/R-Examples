# https://github.com/dmpe/urlshorteneR/issues/9

library(testthat)
library(urlshorteneR)
library(httr)
library(jsonlite)
library(stringr)

bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8",
                          secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")

test_that("issue 9 is fixed", {
  expect_error(expect_message(bitly_LinksShorten(longUrl = ""), "MISSING_ARG_URI"))
})
