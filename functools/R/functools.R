#' functools: Extending Functional Programming in R
#'
#' functools extends functional programming in R. It has three main goals:
#'
#' \itemize{
#' \item Add support to the usual higher order functional suspects (Map, Reduce, Filter, etc.) without extending any core R objects.
#' \item Use a consistent API to access different functionals in base R such as `lapply` or `apply`.
#' \item Provide blazing fast performance for in-memory data by writing key pieces in C++ and options for parallelization, where possible.
#' }
#'
#' functools achieves these goals through three main types of function design patterns:
#'
#' \itemize{
#' \item Closures (functions that take data and return functions)
#' \item Functionals (functions that take functions and return data)
#' \item Function Operators (functions that take functions and return functions)
#' }
#'
#' To learn more about functools, start with the vignettes:
#' \code{browseVignettes(package = "functools")}
#'
#' @docType package
#' @name functools
NULL
