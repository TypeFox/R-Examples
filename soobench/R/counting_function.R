##' Return a new function which is identical to the \code{soofunction}
##' passed in except that all function evaluations are counted.
##'
##' @param fun A test function (class \code{soo_function}).
##'
##' @examples
##' f <- counting_function(double_sum_function(5))
##' number_of_evaluations(f)
##' y <- f(random_parameters(1, f))
##' number_of_evaluations(f)
##' reset_evaluation_counter(f)
##' number_of_evaluations(f)
##' y <- f(random_parameters(21, f))
##' number_of_evaluations(f)
##'
##' @seealso \code{\link{number_of_evaluations}},
##'          \code{\link{reset_evaluation_counter}}
##' @export
counting_function <- function(fun) {            
  force(fun)
  stopifnot("soo_function" %in% class(fun))
  if ("counting_function" %in% class(fun))
    stop("Function already is of type 'counting_function'.")

  count <- 0L
  cfun <- function(x, ...) {
    count <<- count + if (is.matrix(x)) ncol(x) else 1L
    fun(x, ...)
  }
  attributes(cfun) <- attributes(fun)
  class(cfun) <- c("counting_function", class(fun))
  cfun
}

##' Return the number of times a test function has been evaluated.
##'
##' The test function must be wrapped by
##' \code{\link{counting_function}} for this to work.
##' 
##' @param fun A \code{counting_function} as returned by
##'   \code{\link{counting_function}}.
##' 
##' @return The current value of the evaluation counter.
##'
##' @export
number_of_evaluations <- function(fun) {
  stopifnot("counting_function" %in% class(fun))
  environment(fun)$count
}

##' Reset the evaluation counter of a test function.
##'
##' The test function must be wrapped by
##' \code{\link{counting_function}} for this to work.
##' 
##' @param fun A \code{counting_function} as returned by
##'   \code{\link{counting_function}}.
##'
##' @return The current value of the evaluation counter.
##' 
##' @export
reset_evaluation_counter <- function(fun) {
  stopifnot("counting_function" %in% class(fun))
  last_count <- environment(fun)$count
  environment(fun)$count <- 0L
  last_count
}
