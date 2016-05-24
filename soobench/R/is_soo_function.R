##' Check if a function is a SOO function.
##'
##' @param fn Function to check.
##' @return \code{TRUE} if \code{fn} is a proper SOO function, else \code{FALSE}.
##' @export
is_soo_function <- function(fn)
  inherits(fn, "soo_function")
