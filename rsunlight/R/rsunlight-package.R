#' Sunlight Foundation data from R
#'
#' You need API keys for Sunlight Foundation APIs. Please get your own API keys if you
#' plant to use these functions for Sunlight Labs (http://services.sunlightlabs.com/).
#' We set up the functions so that you can use either env vars, or R options. For env
#' vars, put an entry in your \code{.Renviron} file with the name \code{SUNLIGHT_LABS_KEY},
#' so the full thing would be \code{SUNLIGHT_LABS_KEY=<key>}. For R options, put the key in
#' your \code{.Rprofile} file like \code{options(SunlightLabsKey = "key")}. Both are called
#' on R startup, and then you don't have to enter your API key for each run of a function.
#'
#' Currently we have functions to interface with the following Sunlight Foundation APIs,
#' where functions for each API are prefixed with a two letter code indicating the service.
#'
#' \itemize{
#'  \item Congress API (`cg`)
#'  \item Open States API (`os`)
#'  \item Capitol Words API (`cw`)
#' }
#'
#' Note that Puerto Rico is not included in Sunlight Foundation data.
#'
#' @importFrom methods is
#' @importFrom utils head
#' @importFrom httr GET content stop_for_status parse_url
#' @importFrom jsonlite fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom stringr str_sub
#' @name rsunlight-package
#' @aliases rsunlight
#' @docType package
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @author Thomas J. Leeper \email{thosjleeper@@gmail.com}
#' @keywords package
NULL

#' Congressional sessions, up to the 114th (2015-2016)
#' @name sessions
#' @docType data
#' @keywords datasets
NULL
