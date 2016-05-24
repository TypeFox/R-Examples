# :vim set ft=R
#' Call a function until it succeeds
#'
#' Designed for accessing network resources that are unreliable,
#' ntry will call a function up to n times, returning its result
#' or fail.
#'
#' @section Usage:
#' ntry(fn, n) \%::\% Function : numeric : .
#'
#' ntry(fn, n)
#'
#' @section Details:
#' Imagine a function that attempts to access a network resource, like
#' a web service or database. Sometimes there will be network timeouts
#' that you want to recover from, without having to change the application
#' logic. This function allows you to do that, while specifying a maximum
#' retry limit so as not to block the function forever.
#'
#' This higher-order function will call the specified function up to n
#' times, returning on the first successful call. The function fn
#' is a closure that takes a single argument representing the attempt number.
#'
#' If calling the function fails n times, then ntry will fail with an error.
#'
#' @name ntry
#' @param fn A single argument function
#' @param n The number of attempts to call the function
#' @return The result of calling fn
#'
#' @examples
#' \dontrun{
#' fn <- function(i) {
#'   x <- sample(1:4, 1)
#'   flog.info("x = %s",x)
#'   if (x < 4) stop('stop') else x
#' }
#' ntry(fn, 6)
#' }
ntry(fn, n) %::% Function : numeric : .
ntry(fn, n) %as% {
  callCC(function(xfn) {
    lapply(1:n, function(i) tryCatch(xfn(fn(i)), error=function(e) NULL))
    stop(sprintf("ntry failed after %s attempts",n))
  })
}
