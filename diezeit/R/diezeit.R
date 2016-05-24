#' @title ZEIT ONLINE Content API
#' @description A wrapper for the ZEIT ONLINE Content API, 
#' available at \url{http://developer.zeit.de}. It gives access to articles 
#' and corresponding metadata from the ZEIT archive and from ZEIT ONLINE. 
#' A personal API key is required for usage.
#' @name diezeit
#' @docType package
#' @details Accessing the ZEIT archive requires an API key, that can be requested at http://developer.zeit.de/quickstart. Registration is free and allows for API-Access with a limit of 10,000 requests per day. If you do not want to enter your key for each R session, put the following in your .Renviron or .Rprofile file:
#' \code{ZEIT_KEY=PUTYOURKEYHERE}
#' @seealso \code{\link{zeit_client}} for client information and usage,
#' \code{\link{zeit_search}} for ZEIT archive search or \code{\link{zeit_get}} 
#' to get content from the ZEIT archive.
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom jsonlite fromJSON
#' @import brew
#' @import grDevices
#' @import methods
#' @import utils
#' @aliases diezeit diezeit-package
NULL
