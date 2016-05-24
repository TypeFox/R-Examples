## ---- echo = FALSE-------------------------------------------------------
library(httr)
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## ------------------------------------------------------------------------
github_GET <- function(path, ..., pat = github_pat()) {
  auth <- github_auth(pat)
  req <- GET("https://api.github.com", path = path, auth, ...)
  github_check(req)

  req
}

github_POST <- function(path, body, ..., pat = github_pat()) {
  auth <- github_auth(pat)

  stopifnot(is.list(body))
  body_json <- jsonlite::toJSON(body)

  req <- POST("https://api.github.com", path = path, body = body_json,
    auth, post, ...)
  github_check(req)

  req
}

## ------------------------------------------------------------------------
github_auth <- function(pat = github_pat()) {
  authenticate(pat, "x-oauth-basic", "basic")
}

github_check <- function(req) {
  if (req$status_code < 400) return(invisible())

  message <- github_parse(req)$message
  stop("HTTP failure: ", req$status_code, "\n", message, call. = FALSE)
}

github_parse <- function(req) {
  text <- content(req, as = "text")
  if (identical(text, "")) stop("No output to parse", call. = FALSE)
  jsonlite::fromJSON(text, simplifyVector = FALSE)
}

github_pat <- function() {
  Sys.getenv('GITHUB_PAT')
}

has_pat <- function() !identical(github_pat(), "")

## ------------------------------------------------------------------------
rate_limit <- function() {
  req <- github_GET("rate_limit")
  github_parse(req)
}

if (has_pat()) {
  str(rate_limit())
}

## ------------------------------------------------------------------------
rate_limit <- function() {
  req <- github_GET("rate_limit")
  core <- github_parse(req)$resources$core

  reset <- as.POSIXct(core$reset, origin = "1970-01-01")
  cat(core$remaining, " / ", core$limit,
    " (Reset ", strftime(reset, "%H:%M:%S"), ")\n", sep = "")
}

if (has_pat()) {
  rate_limit()
}

## ------------------------------------------------------------------------
github_parse <- function(req) {
  text <- content(req, as = "text")
  if (identical(text, "")) stop("")
  jsonlite::fromJSON(text, simplifyVector = FALSE)
}

## ------------------------------------------------------------------------
github_parse <- function(req) {
  text <- content(req, as = "text")
  if (identical(text, "")) stop("No output to parse", call. = FALSE)
  jsonlite::fromJSON(text, simplifyVector = FALSE)
}

## ------------------------------------------------------------------------
authenticate("username", "password")

## ------------------------------------------------------------------------
authenticate("ddfa3d40d5855d6ba76b7003fd4", "")

## ------------------------------------------------------------------------
github_pat <- function(force = FALSE) {
  env <- Sys.getenv('GITHUB_PAT')
  if (!identical(env, "") && !force) return(env)

  if (!interactive()) {
    stop("Please set env var GITHUB_PAT to your github personal access token",
      call. = FALSE)
  }

  message("Couldn't find env var GITHUB_PAT. See ?github_pat for more details.")
  message("Please enter your PAT and press enter:")
  pat <- readline(": ")

  if (identical(pat, "")) {
    stop("Github personal access token entry failed", call. = FALSE)
  }

  message("Updating GITHUB_PAT env var to PAT")
  Sys.setenv(GITHUB_PAT = pat)

  pat
}

## ------------------------------------------------------------------------
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(purl = NOT_CRAN)

