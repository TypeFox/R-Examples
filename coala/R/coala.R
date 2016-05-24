#' A Framework for Coalescent Simulation in R
#'
#' This package allows to specify and simulate coalescent models from
#' within R. The `introduction` vignette is a good place to start.
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib coala
#' @importFrom assertthat assert_that
"_PACKAGE"

# Mute warnings about R6 object internals
#' @importFrom utils suppressForeignCheck
suppressForeignCheck(c("self", "private", "super"))

release_questions <- function() {
  c("Have you tested the package with valgrind?")
}
