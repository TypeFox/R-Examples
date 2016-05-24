# the LongURL api requires the use of a custom user agent
LONGURL_USER_AGENT <- "longurl-rstats-pkg"

# timeout
LONGURL_TIMEOUT <- 1.5

# this is the base endpoint for the LongURL API
LONGURL_ENDPOINT <- "http://api.longurl.org/v2/%s"

#' Retrieve all the URL shortener services known to the 'LongURL' API
#'
#' @export
#' @return \code{data_frame} (compatible with \code{data.frame}) of results
#'        with the the short URL TLD in \code{tld} and a regular expression
#'        of compatible URLs in \code{regex}). Will populate the returned
#'        \code{data.frame} with \code{NA} if there are severe issues connecting
#'         to the LongURL API.
#' @examples
#' short_svcs <- known_services()
#' head(short_svcs)
known_services <- function() {

  url <- sprintf(LONGURL_ENDPOINT, "services")

  resp <- tryCatch(GET(url, query=list(format="json"),
                       user_agent(LONGURL_USER_AGENT),
                       timeout(LONGURL_TIMEOUT)),
                   error=function(err) return(NA_character_))

  if (is.na(resp)) return(data.frame(domain=NA_character_, regex=NA_character_))

  warn_for_status(resp)

  tmp <- content(resp)

  data_frame(domain=as.vector(sapply(tmp, "[[", "domain", USE.NAMES=FALSE)),
             regex=as.vector(sapply(tmp, "[[", "regex", USE.NAMES=FALSE)))

}

#' Expand a vector of (short) URLs using the longurl service
#'
#' Pass in a vector of URLs (ostensibly "short" URLs) and receive
#' a \code{data_frame} of the original URLs and expanded URLs via the
#' 'LongURL' service.
#'
#' @param urls_to_expand character vector of URLs
#' @param check run an extra \code{HEAD} request on the expanded URL to determine
#'        validity. This is an expensive operation, so recommended usage is to run
#'        this only on URLs that did not seem to expand.
#' @param warn show any warnings (API or otherwise) as messages
#' @param .progress display a progress bar (generally only useful in
#'        interactive sesions)
#' @return \code{data_frame} (compatible with \code{data.frame}) of results
#'        with the orignial URLs in \code{orig_url} and expanded URLs in
#'        \code{expanded_url}). Will return \code{NA} if
#'        there are severe issues connecting to the LongURL API.
#' @export
#' @examples
#' test_urls <- c("http://t.co/D4C7aWYIiA",
#'                "ift.tt/1L2Llfr")
#' big_urls <- expand_urls(test_urls)
#' head(big_urls)
expand_urls <- function(urls_to_expand, check=FALSE, warn=TRUE,
                        .progress=interactive()) {

  doapply <- ifelse(.progress, pbsapply, sapply)

  data_frame(orig_url=urls_to_expand,
             expanded_url=doapply(urls_to_expand, expand_url,
                                  check=check, warn=warn, USE.NAMES=FALSE))

}

#' the thing that does all the work
#' @noRd
expand_url <- function(url_to_expand, check=FALSE, warn=TRUE) {

  # make the API URL
  url <- sprintf(LONGURL_ENDPOINT, "expand")

  # use the API
  resp <- tryCatch(GET(url, query=list(url=url_to_expand,
                                       format="json"),
                       user_agent(LONGURL_USER_AGENT),
                       timeout(LONGURL_TIMEOUT)),
                   error=function(err) { return(NA_character_) })

  if (is.na(resp)) return(NA_character_)

  # warn for API errors
  if (warn) warn_for_status(resp)

  # response object
  tmp <- content(resp)

  # if bad response and/or API long-url not populated kick back NA
  if ((resp$status != 200) | (!(names(tmp) %in% c("long-url")))) return(NA)

  # _expensive_ validity check of expanded URL
  if (check) {
    chk <- HEAD(url_to_expand)
    if (warn) warn_for_status(chk)
    if (chk$status != 200) return(NA)
  }

  return(tmp$`long-url`)

}
