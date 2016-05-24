## ------------------------------------------------------------------------
library(urlshorteneR)
# print(getwd())
# 
# btoken_path <- file.path("..", "tests", "testthat", "bitly_token.rds")
# gtoken_path <- file.path("..", "tests", "testthat", "googl_token.rds")
 
if (interactive()) {

bitly_token <- readRDS("../tests/testthat/bitly_token.rds")
googl_token <- readRDS("../tests/testthat/googl_token.rds")

# You should register a new pair of keys yourself. DO NOT USE MINE as this may not work. 
# bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8", secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")

bitly_UserMetricsPopularLinks(unit = "month", units = -1, limit = 100)
}

## ------------------------------------------------------------------------
if (interactive()) {

bitly_LinksMetricsEncodersByCount(link = "http://bit.ly/DPetrov", my_network = "false", limit = 100)
}

## ------------------------------------------------------------------------
if (interactive()) {

bitly_UserInfo()

bitly_UserTrackingDomains()
}

## ------------------------------------------------------------------------
if (interactive()) {

bitly_IsProDomain(domain = "nyti.ms")
}

## ------------------------------------------------------------------------
if (interactive()) {

# googl_auth(key = "806673580943-78jdskus76fu7r0m21erihqtltcka29i.apps.googleusercontent.com", secret = "qItL-PZnm8GFxUOYM0zPVr_t")

g2 <- googl_LinksShorten(longUrl = "https://developers.google.com/url-shortener/v1/url/insert")
g2
g1 <- googl_LinksExpand(shortUrl = "http://goo.gl/vM0w4", showRequestURL = F)
g1
}

## ------------------------------------------------------------------------
isgd_LinksShorten(longUrl = "http://debil.cz/", showRequestURL = TRUE)

