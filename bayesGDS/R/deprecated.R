## deprecated.R -- Part of the bayesGDS package
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.



## Place to hold functions that will not be maintained in the future

#' @name Deprecated
#' @aliases vech inv.vech
#' @title Deprecated functions
#' @description These functions were in earlier versions, but will no
#' longer be maintained in this package.  They will likely be moved to
#' another package a some time.
NULL


#' @title vech operator on a square matrix
#' @param M a matrix
#' @return A vector containing the lower triangle of M, ordered column-wise.
#' @rdname Deprecated
#' @export
vech <- function( M )
{
    .Deprecated("matrix")

    if ( nrow(M)!=ncol(M))
        stop( "argument M is not a square numeric matrix" )
    return( t( t( M[!upper.tri(M)] ) ) )
}

#' @title inverse vech operator on a vector
#' @param y A vector of conforming length
#' @return A k x k lower triangular matrix
#' @rdname Deprecated
#' @export
inv.vech <- function( y ) {

    .Deprecated("as.vector")

    n <- length(y)
    R <- (sqrt(1+8*n)-1)/2

    if (R != as.integer(R)) {
        stop ("in function inv.vech:  vector will not fit in square matrix\n")
    }

    res <- matrix(0,R, R)
    res[!upper.tri(res)] <- y
    return(res)
}


#' @title Logit transformation
#' @param p A scalar, vector or matrix, where each element is between
#' 0 and 1.
#' @return result = log(p/(1-p))
#' @rdname Deprecated
#' @export
logit <- function(p) {

  if (any(p<=0) || any(p>=1)) {
    stop (" in function logit:  all elements must be in open interval (0,1)")
  }

  res <- log(p) - log(1-p)
  return(res)

}

#' @title Inverse logit transformation
#' @param x A scalar, vector or matrix
#' @return result = exp(x)/(1+exp(x))
#' @rdname Deprecated
#' @export
inv.logit <- function(x) {

  ## numerically stable inv.logit

  w.max <- x>=log(.Machine$double.xmax)

  res <- exp(x - log1p(exp(x)))
  res[w.max] <- 1
  return(res)

}

#' @title Log inverse logit transformation
#' @return result = log[exp(x)/(1+exp(x))]
#' @rdname Deprecated
#' @export
log_inv.logit <- function(x) {
    w.max <- x>=log(.Machine$double.xmax)
    w.min <- x<=log(.Machine$double.xmin)
    ww <- !(w.min | w.max)

    res <- x
    res[ww] <- x[ww] - log1p(exp(x[ww]))
    res[w.max] <- 0
    return(res)
}

