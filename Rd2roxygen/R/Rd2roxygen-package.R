#' Convert Rd to roxygen documentation and utilities to enhance documentation
#' and building packages
#'
#' This package contains functions to convert Rd to roxygen documentation. It
#' can parse an Rd file to a list (\code{\link{parse_file}}), create the roxygen
#' documentation (\code{\link{create_roxygen}}) and update the original R script
#' (e.g. the one containing the definition of the function) accordingly
#' (\code{\link{Rd2roxygen}}). This package also provides utilities which can
#' help developers build packages using roxygen more easily (\code{\link{rab}}).
#'
#' @name Rd2roxygen-package
#' @docType package
#' @author Hadley Wickham and Yihui Xie
#' @note There is no guarantee to generate perfect roxygen comments that can be
#'   converted back to the original Rd files. Usually manual manipulations on
#'   the roxygen comments are required. For example, currently `@@S3method` is
#'   not included in the comments, and `@@rdname` is not supported either (users
#'   have to move the roxygen comments around and add the appropriate tags by
#'   themselves). Patches (as pull requests) are welcome through GitHub:
#'   \url{https://github.com/yihui/Rd2roxygen/}.
#'
#'   This package is not thoroughly tested, so it is likely that it fails to
#'   convert certain parts of Rd files to roxygen comments. As mentioned before,
#'   you have to manually deal with these problems. You are welcome to report
#'   other serious issues via \url{https://github.com/yihui/Rd2roxygen/issues}.
#' @importFrom formatR tidy_source
#' @examples
#' ## see the package vignette: vignette('Rd2roxygen')
NULL
