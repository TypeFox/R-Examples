#' Rfit Internal Print Functions
#' 
#' These functions print the output in a user-friendly manner using the
#' internal R function \code{print}.
#' 
#' 
#' @aliases print.rfit print.summary.rfit print.drop.test print.oneway.rfit
#' print.summary.oneway.rfit print.raov
#' @param x An object to be printed
#' @param digits number of digits to display
#' @param \dots additional arguments to be passed to \code{print}
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\link{rfit}}, \code{\link{summary.rfit}},
#' \code{\link{drop.test}}
#' @export print.rfit
print.rfit <- function (x, ...) {
    cat("Call:\n")
    print(x$call)
    coef <- coef(x)
    cat("\nCoefficients:\n")
    print(coef, ...)
}
