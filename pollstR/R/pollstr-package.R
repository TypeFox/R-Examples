#' R client for the Huffpost Pollster API
#'
#' This pacakge provides an R interface to the Huffington Post Pollster API.
#' Pollster provides programmatic access to opinion polls collected by the Huffington Post.
#'
#' See \url{http://elections.huffingtonpost.com/pollster/api} for more details on the API.
#'
#' @name pollstR
#' @docType package
#' @import httr
#' @import plyr
#' @import jsonlite
#' @importFrom utils head
NULL

.POLLSTR_API_URL <- "http://elections.huffingtonpost.com/pollster/api"

get_url <- function(url, as = "parsed") {
    response <- GET(url)
    stop_for_status(response)
    content(response, as = as)
}

convert_df <- function(x) {
    for (i in names(x)) {
        if (is.null(x[[i]])) {
            x[[i]] <- NA
        }
    }
    data.frame(x, stringsAsFactors = FALSE)
}

# election date entry
electiondate2date <- function(x) {
    if (is.null(x)) {
        as.Date(NA_character_)
    } else {
        as.Date(x, "%Y-%m-%d")
    }
}
