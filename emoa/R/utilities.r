##
## utilities.r - Internal utility functions
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

##' @useDynLib emoa do_which_points_on_edge
NA

rsbxbeta <- function(n, nc) {
    ## U ~ [0, 1] , twou := 2 * U ~ [0, 2]
    twou <- runif(n, 0, 2.0)

    e <- 1/(nc + 1)
    beta <- ifelse(twou < 1, twou, 1/(2 - twou))^e
    return(beta)
}

##' Clip value to a given range
##'
##' Clip \eqn{x} to the interval \eqn{[l, u]}. This is useful to enforce
##' box constraints.
##'
##' @param x Value to clip.
##' @param l Lower limit.
##' @param u Upper limit.
##'
##' @return l if x < l, u if x > u else x.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @export
inbounds <- function(x, l, u) {
  ifelse(x < l, l, ifelse(x > u, u, x))
}

##' Return first non null argument.
##'
##' This function is useful when processing complex arguments with multiple
##' possible defaults based on other arguments that may or may not have been
##' provided.
##'
##' @param ... List of values.
##' @return First non null element in \code{...}.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @export
coalesce <- function(...) {
  l <- list(...)
  isnull <- sapply(l, is.null)
  l[[which.min(isnull)]]
}
