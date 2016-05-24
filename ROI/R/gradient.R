## gradients
## code based on a patch submitted by Olaf Mersmann.
## slightly modified

##' Extract the gradient from its argument (typically a ROI
##' object of class \code{"objective"}).
##'
##' @title Extract Gradient information
##' @param x an object used to select the method.
##' @param \ldots further arguments passed down to the
##' \code{\link[numDeriv]{grad}()} function for calculating gradients (only for \code{"F_objective"}).
##' @return a \code{"function"}.
##' @author Stefan Theussl
##' @export
G <- function( x, ... )
    UseMethod("G")

G.F_objective <- function( x, ... ){
    args <- list(...)
    args$func <- x$F
    g <- x$G
    if(is.null(g))
        g <- function(x){
            args$x <- x
            do.call(numDeriv::grad, args = args)
        }
    stopifnot( is.function(g) )
    g
}

## FIXME: as.numeric method for stms?
G.L_objective <- function( x, ... ){
    L <- terms(x)$L
    function(x)
        as.numeric(as.matrix(L))
}

G.Q_objective <- function( x, ... ){
    L <- terms(x)$L
    Q <- terms(x)$Q
    function(x)
        c(slam::tcrossprod_simple_triplet_matrix(Q, t(x)) + L)
}

#G.function <- function(x)

#G.default <- function(x)
#    NULL
