library(testthat)
library(urlshorteneR)
library(httr)
library(jsonlite)
library(stringr)

googl_token <- googl_auth(key = "806673580943-78jdskus76fu7r0m21erihqtltcka29i.apps.googleusercontent.com",
             secret = "qItL-PZnm8GFxUOYM0zPVr_t")
bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8",
             secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")

# bitly_token <- readRDS("../bitly_token.rds")

context("User Metrics")

test_that("Returns aggregate metrics about the countries referring click traffic to all of the authenticated user's Bitlinks.", {
  umcoun <- bitly_UserMetricsCountries(unit = "day", units = -1, limit = 100, rollup = "true")
  expect_named(umcoun, c("country","clicks")) 
})

test_that("Returns the aggregate number of clicks on all of the authenticated user's Bitlinks.", {
  umc <- bitly_UserMetricsClicks(unit = "day", units = -1, limit = 100, rollup = "true")
  expect_more_than(umc, 5)
  umcc <- bitly_UserMetricsClicks(unit = "day", units = -1, limit = 100, rollup = "false")
  expect_named(umcc, c("dt", "clicks")) 
})

test_that("Returns the authenticated user's most-clicked Bitlinks (ordered by number of clicks) in a given time period.", {
  umpl <- bitly_UserMetricsPopularLinks(unit = "month", units = -1, limit = 100)
  expect_named(umpl, c("link", "clicks")) 
})

test_that("Returns aggregate metrics about the pages referring click traffic to all of the authenticated user's Bitlinks.", {
  umrr <- bitly_UserMetricsReferrers(unit = "day", units = -1, limit = 100, rollup = "true")
  expect_named(umrr, c("referrer", "clicks")) 
})

test_that("Returns aggregate metrics about the domains referring click traffic to all of the authenticated user's Bitlinks.", {
  umrd <- bitly_UserMetricsReferringDomains(unit = "day", units = -1, limit = 100, rollup = "true", exclude_social_networks = "false")
  expect_named(umrd, c("domain", "clicks"))
  expect_message(bitly_UserMetricsReferringDomains(unit = "day", units = -1, limit = 100, rollup = "true", exclude_social_networks = "true")
                 , "You have zero referring domains given your function input.") 
})

test_that("Returns the number of Bitlinks created in a given time period by the authenticated user.", {
  umsc <- bitly_UserMetricsShortenCounts(unit = "day", units = -1, limit = 100, rollup = "true")
  expect_more_than(umsc, 5)
  umscf <- bitly_UserMetricsShortenCounts(unit = "day", units = -1, limit = 100, rollup = "false")
  expect_named(umscf, c("dt", "shortens"))
})

