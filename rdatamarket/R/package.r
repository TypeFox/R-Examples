#' Data access API for DataMarket.com
#'
#' Interacts with DataMarket.com, fetching data in timeseries (\code{zoo}) and
#' long form (suitable for \code{reshape}), and metadata on datasets.
#'
#' @section Quick start:
#'
#' Just find the data you want on \url{http://datamarket.com/} and then
#' paste the URL from your browser into \code{\link{dmlist}} or
#' \code{\link{dmseries}}:
#'
#' \code{plot(dmseries("http://datamarket.com/data/set/zy0/#ds=zy0!8o7=2"))}
#'
#' \code{plot(dmseries("http://data.is/nyFeP9"))}
#'
#' \code{l <- dmlist("http://data.is/nyFeP9")}
#'
#' or into \code{\link{dminfo}} to get the dataset metadata:
#'
#' \code{ds <- dminfo("http://data.is/nyFeP9")}
#'
#' @section DS parameters:
#'
#' In functions that take a parameter named \code{ds}, that parameter can be:
#' \itemize{
#'   \item a dataset ID (\code{"17tm"})
#'   \item a DS string (\code{"17tm!kqc=a"}) consisting of a dataset ID
#'         along with some dimension filtering specifications; for more about
#'         DS strings, see Datamarket API documentation at
#'         \url{http://datamarket.com/api/v1/#about-ds}
#'   \item a query string (\code{"foo=bar&ds=17tm!kqc=a&baz=xyzzy"})
#'         containing such a string in a \code{ds} parameter
#'   \item a URL (\code{"http://datamarket.com/data/set/17tm/#ds=17tm!kqc=a"})
#'         containing such a query string after \code{?} or \code{#}
#'   \item such a URL represented by a redirect (\code{"http://data.is/r6JEsC"})
#'         in a short-URL service (\code{data.is, bit.ly, is.gd, t.co, url.is})
#'   \item a dataset object, obtained by calling \code{\link{dminfo}}.
#' }
#'
#' Thus when viewing any dataset on datamarket.com, you can copy the URL from
#' your browser bar and paste into a call to \code{\link{dmseries}} or
#' \code{\link{dmlist}} to load the same data into R, or a call to
#' \code{\link{dminfo}} to get a dataset object with which to make data queries.
#'
#' @section Dimension filtering:
#'
#' The functions \code{\link{dmseries}} and \code{\link{dmlist}} can accept
#' parameters named for the dimensions of the dataset, in order to filter the
#' data: \code{dmlist("17tm", Country="Algeria")}
#'
#' @import zoo RCurl RJSONIO
#' @docType package
#' @name rdatamarket-package
#' @example demo/rdatamarket-usage.r
NULL

.rdatamarketEnv <- new.env()
.rdatamarketEnv$curlopts <- curlOptions()
.rdatamarketEnv$api.key <- NULL
