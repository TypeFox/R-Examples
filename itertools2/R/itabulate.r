#' Iterator that maps a function to a sequence of numeric values
#'
#' Constructs an iterator that maps a given function over an indefinite sequence
#' of numeric values. The input the function \code{f} is expected to accept a
#' single numeric argument. The sequence of arguments passed to \code{f} begin
#' with \code{start} and are incremented by \code{step}.
#'
#' @export
#' @param f the function to apply
#' @param start sequence's initial value
#' @param step sequence's step size
#' @return an iterator that returns the mapped values from the sequence
#'
#' @examples
#' it <- itabulate(f=function(x) x + 1)
#' take(it, 4) # 2 3 4 5
#' 
#' it2 <- itabulate(f=function(x) x^2, start=-3)
#' take(it2, 6) # 9 4 1 0 1 4
#'
#' it3 <- itabulate(abs, start=-5, step=2)
#' take(it3, 6) # 5 3 1 1 3 5
#'
#' it4 <- itabulate(exp, start=6, step=-2)
#' take(it4, 4) # exp(c(6, 4, 2, 0))
#'
itabulate <- function(f, start=1, step=1) {
  start <- as.numeric(start)
  step <- as.numeric(step)
  
  if (length(start) != 1) {
    stop("'start' must be a numeric value of length 1")
  }
  if (length(step) != 1) {
    stop("'step' must be a numeric value of length 1")
  }

  imap(f=f, icount(start=start, step=step))
}
