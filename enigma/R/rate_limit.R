#' Get rate limit data.
#'
#' @export
#'
#' @param key (character) Required. An Enigma API key. Supply in the function call, or store in
#' your \code{.Rprofile} file, or do \code{options(enigmaKey = "<your key>")}. Obtain an API key
#' by creating an account with Enigma at \url{http://enigma.io}, then obtain an API key from
#' your account page.
#' @param ... Named options passed on to \code{\link[httr]{GET}}
#' @examples \dontrun{
#' rate_limit()
#' }

rate_limit <- function(key=NULL, ...) {
  key <- check_key(key)
  url <- sprintf('%s/limits/%s', en_base(), key)
  json <- enigma_GET(url, NULL, ...)
  structure(json, class = "enigma_rate")
}

#' @export
print.enigma_rate <- function(x, ...) {
  cat("Period: ", x$period, "\n", sep = "")
  cat("Rate limit: 5000", "\n", sep = "")
  cat("Remaining - data:  ", x$data, "\n", sep = "")
  cat("Remaining - meta:  ", x$meta, "\n", sep = "")
  cat("Remaining - export:  ", x$export, "\n", sep = "")
  cat("Resets in:  ", time(x$seconds_remaining), "\n", sep = "")
}

time <- function(x) {
  x <- as.integer(x)
  if (x > 3600) {
    paste0(x %/% 3600, " hours")
  }
  else if (x > 300) {
    paste0(x %/% 60, " minutes")
  }
  else if (x > 60) {
    paste0(round(x/60, 1), " minutes")
  }
  else {
    paste0(x, "s")
  }
}
