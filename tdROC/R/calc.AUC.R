#' Calculate the area under a ROC curve (AUC)
#'
#' This function reads in a vector of sensitivity and a vector of specificity and calculates
#' the area under the curve (AUC) by trapezoidal integration.
#'@param sens a numerical vector of sensitivity values within the range of (0, 1).
#'@param spec a numerical vector of specificity values within the range of (0, 1).
#'@return AUC as a numerical scalar.
#'@note This function sorts \code{sens} and \code{1-spec} in an increasing order.
#'      A 0 and 1 will be added to the two ends of the sorted vectors. The Area Under the Curve (AUC) is obtained by trapezoidal
#'       integration of the area under the piecewise linear curve obtained by connecting
#'        points in \code{sens} and \code{1-spec}. .
#'@examples
#'sens <- c(0.99, 0.97, 0.83, 0.60, 0.40, 0.20) ;
#'spec <- c(0.50, 0.61, 0.80, 0.90, 0.95, 0.98) ;
#'plot( 1-spec, sens, type = "l" ) ;
#'points(1-spec, sens) ;
#'calc.AUC( sens, spec ) ;
#'
#'@export


calc.AUC <- function( sens, spec ) {
  # Given a vector of sensitivity and specificity, calculate the area under the
  # ROC curve (AUC)
  # Arguments:
  #  -- sens: sensitivity
  #  -- spec: specificity
  # Return:
  #  -- AUC as a numerical scalar
  # NOTE: The AUC is obtained by trapezoidal integration of the area under the
  #       piecewise linear curve obtained by connecting the sensitivity and
  #       specificity.
  o <- order(sens, 1-spec) ;
  y <- sens[o] ;
  x <- 1-spec[o] ;
  x <- c(0, x, 1) ;
  y <- c(0, y, 1) ;
  m <- length(x) ;
  x1 <- x[-m] ;
  x2 <- x[-1] ;
  y1 <- y[-m] ;
  y2 <- y[-1] ;
  AUC <- sum( (y1 + y2)*(x2 - x1)/2 ) ;
  AUC ;
}

