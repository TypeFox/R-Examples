##' Wrap vector
##'
##' Wrap vector
##' @param x Vector or integer 
##' @param delta Shift
##' @param ... Additional parameters
##' @export
##' @examples
##' wrapvec(5,2)
wrapvec <- function(x,delta=0L,...) {
    if (length(x)==1 && floor(x)==x && x>0) {
        x <- seq(x)
    }
    if (delta==0L) return(x)
    x[(seq_along(x)+delta-1L)%%length(x)+1L]
}
