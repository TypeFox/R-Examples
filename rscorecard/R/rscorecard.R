#' rscorecard: A Method to Download College Scorecard Data.
#'
#' The rscorecard package provides a series of piped functions (a la
#' \href{http://cran.r-project.org/package=dplyr}{dplyr}) to
#' facilitate downloading Department of Education College Scorecard
#' data.  In reality it is simply a method for converting idiomatic R
#' code into a properly formatted URL string that is then
#' queried. This package requires an API key, which can be requested
#' at \url{https://api.data.gov/signup/}.
#'
#' All command pipes must start with \code{sc_init()}, end with
#' \code{sc_get()}, and be linked with the magrittr pipe function,
#' \code{\%>\%}.  Internal commands, \code{\link{sc_select}},
#' \code{\link{sc_filter}}, \code{\link{sc_year}},
#' \code{\link{sc_zip}}, can come in any order in the pipe chain. Only
#' \code{\link{sc_select}} is required.
#'
#' @docType package
#' @name rscorecard
NULL
