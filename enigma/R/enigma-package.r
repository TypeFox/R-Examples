#' enigma, an R client for Enigma.io
#'
#' Enigma holds government data and provides a really nice set of APIs for data, metadata, and 
#' stats on each of the datasets. That is, you can request a dataset itself, metadata on the 
#' dataset, and summary statistics on the columns of each dataset.
#'
#' The functions:
#'
#' \itemize{
#'  \item \code{\link[enigma]{enigma_data}} - Fetch and dataset, and filter on columns or rows.
#'  \item \code{\link[enigma]{enigma_metadata}} - Get metadata on datasets.
#'  \item \code{\link[enigma]{enigma_stats}} - Get columnwise statistics on datasets.
#'  \item \code{\link[enigma]{enigma_fetch}} - Get gzipped csv of a dataset. Goes along with 
#'  \code{\link[enigma]{enigma_read}}
#'  \item \code{\link[enigma]{rate_limit}} - Get columnwise statistics on datasets.
#' }
#' 
#' An API key is required to use this package. You can supply your key in each function call, or 
#' store in your key in your \code{.Rprofile} file, or execute 
#' \code{options(enigmaKey = "<your key>")} in your R console. Obtain an API key 
#' by creating an account with Enigma at \url{http://enigma.io}, then get an API key from 
#' your Enigma account page.
#'
#' @importFrom methods is
#' @importFrom utils read.delim head
#' @name enigma
#' @docType package
NULL
