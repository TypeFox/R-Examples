#' Check monotonicity
#'
#' This internel function determines if a numeric vector X is monotone (increasing or decreasing) or not. Ties in
#' the sequence are considered as monotone.
#' @param X a numerical vector to be tested.
#' @return Returns True or False.
#' @examples
#'
#' a <- 1:10 ;
#' b <- c( 1:5, 4:2 ) ;
#' c <- c( 1:4, 4, 4, 4, 5:10 ) ;
#' d <- c( 10:6, 3:-2 ) ;
#'
#' is.monotone( a )
#' is.monotone( b )
#' is.monotone( c )
#' is.monotone( d )
#'
#' @export
is.monotone <- function(X) {
  # determine if a vector x is monotone (increasing or decreasing) or not
  m <- length(X) ;
  if ( all( X[-m] - X[-1] >= 0 ) | all( X[-m] - X[-1] <= 0 ) ) {
    ans <- TRUE ;
  } else {
    ans <- FALSE ;
  }
  ans ;
}
