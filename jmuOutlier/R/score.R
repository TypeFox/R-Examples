score <-
function(x, y=NULL, expon=FALSE) {
  # Generates van der Waerden scores (i.e., normal quantiles) or exponential 
  #    (similar to Savage) scores, for combined data \code{x} and \code{y}.
  # 'x': A positive integer equal to the number of desired scores when \code{y} is \code{NULL},
  #      or \code{x} is a vector of observations.
  # 'y': An optional vector of observations, typically used with two-sample tests.
  # 'expon': Logical; if \code{FALSE} (default), van der Waerden scores are computed, even for ties.
  #                   If \code{TRUE}, Exponential scores are computed, and interpolation is used for ties.
  # Output: Scored values for \code{x}, when \code{y} is \code{NULL}.
  # Output: \code{$x} are the scored values for \code{x}, when \code{y} is not \code{NULL}.
  # Output: \code{$y} are the scored values for \code{y}, when \code{y} is not \code{NULL}.
  # Examples:   score( 10 )
  #             score( 15, expon=TRUE )
  #             score( c(4,7,6,22,13), c(15,16,7) )  # Two samples, including a tie.
  if ( !is.numeric(x) )  stop("'x' must be numeric.")
  if ( !is.numeric(y) & !is.null(y) )  stop("'y' must be numeric or NULL.")
  if ( !is.logical(expon) )  stop("'expon' must be logical.")
  if (is.null(y) & length(x)==1) { if (floor(x)==ceiling(x) & x>=1) {x <- 1:x} }
  n <- length(c(x,y)) ;    if (!expon) {z <- qnorm(rank(c(x,y))/(n+1))}
  if (expon) { u <- cumsum( 1/(n+1-1:n) ) ;   z <- u[ rank(c(x,y)) ]
    for (i in 1:n) { rank0 <- rank(c(x,y))[i] ;   remainder <- rank0 - floor(rank0)
      if (remainder>0) {z[i] <- u[floor(rank0)]*(1-remainder)+u[ceiling(rank0)]*remainder}}}
  if (is.null(y)) { return(z) }
  if (!is.null(y)) { return(structure(list(x=z[1:length(x)], y=z[(length(x)+1):n]))) }
}
