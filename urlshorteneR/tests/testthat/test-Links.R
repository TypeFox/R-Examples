library(testthat)
library(urlshorteneR)
library(httr)
library(jsonlite)
library(stringr)

# bitly_token <- readRDS("bitly_token.rds")

context("Links Bit.ly")

test_that("Query for a Bitlink based on a long URL.", {
  ll <- bitly_LinksLookup(url = "http://www.google.com/")
  expect_equal(dim(ll)[[2]], 2)
})

test_that("Used to return the page title for a given Bitlink.", {
  li <- bitly_LinksInfo(hashIN = "DPetrov", expand_user = "true")
  li2 <- bitly_LinksInfo(shortUrl = "http://bit.ly/DPetrov", expand_user = "false")
  expect_equal(dim(li)[[2]], 11)
  expect_equal(dim(li2)[[2]], 7)
})

test_that("Given a bitly URL or hash (or multiple), returns the target (long) URL.", {
  le <- bitly_LinksExpand(hashIN = "DPetrov")
  le2 <- bitly_LinksExpand(shortUrl = "http://bit.ly/DPetrov")
  expect_equal(dim(le)[[2]], 4)
  expect_named(le2, c("short_url", "long_url", "user_hash", "global_hash"))
})

test_that("Given a bitly URL or hash (or multiple), returns the target (long) URL.", {
  ls <- bitly_LinksShorten(longUrl = "http://slovnik.seznam.cz/")
  ls2 <- bitly_LinksShorten(longUrl = "https://travis-ci.org/dmpe/rbitly/builds/68231423", domain = "j.mp")
  expect_equal(dim(ls)[[2]], 5)
  expect_named(ls2, c("long_url", "url", "hash", "global_hash", "new_hash"))
})

context("Links Goo.gl")

test_that("expanding does work", {
  g1 <- googl_LinksExpand(shortUrl = "http://goo.gl/vM0w4", showRequestURL = F)
  expect_output(g1$original_data$longUrl, "http://www.bi-verdict.com/fileadmin/FreeAnalyses/consolidations.htm")
})

test_that("shorting does work", {
  g2 <- googl_LinksShorten(longUrl = "https://developers.google.com/url-shortener/v1/url/insert", showRequestURL = F)
  expect_more_than(nchar(g2$id), 19)
}) 

