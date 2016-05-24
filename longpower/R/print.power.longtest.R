#' Print method for longitudinal data power calculation object
#' 
#' Print object of class \code{"power.longtest"} in nice layout.
#' 
#' A \code{power.longtest} object is just a named list of numbers and character
#' strings, supplemented with \code{method} and \code{note} elements.  The
#' \code{method} is displayed as a title, the \code{note} as a footnote, and
#' the remaining elements are given in an aligned \sQuote{name = value} format.
#' 
#' @param x Object of class \code{"power.longtest"}.
#' @param \dots further arguments to be passed to or from methods.
#' @return none
#' @seealso \code{\link{liu.liang.linear.power}},
#' \code{\link{diggle.linear.power}}, \code{\link{lmmpower}},
#' @keywords longtest
print.power.longtest <- function(x, ...)
{
  cat("\n    ", x$method, "\n\n")
  note <- x$note
  R <- x$R
  x[c("method","note","R")] <- NULL
  cat(paste(format(names(x), width= 15, justify = "right"),
  format(x), sep= " = "), sep= "\n")
  if(!is.null(note)){
    cat("\n", "NOTE:", note, "\n")
  }else{
    cat("\n", "NOTE: n is the number in *each* group\n")
  }
  if(!is.null(R)){
    cat("\n", "R:\n")
    print(R)
    cat("\n")
  }
  invisible(x)
}
