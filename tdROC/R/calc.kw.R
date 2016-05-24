#' Calculate Kernel weights
#'
#' This function calculate the nearest neighbor kernel weights using uniform kernel weights.
#'
#' @param span a numeric value of the proportion of neighbour observations used, default is 0.1.
#' @param h a numeric value of the bandwidth of kernel weights, defualt is NULL. If not specified, the function used the value of \code{span} to calculate weights.
#'        If both \code{span} and \code{h} are specified, the function will ignore the span and used bandwidth of kernel function instead.
#' @param type a character value of the type of kernel function used to calculate kernel weights. Default is "uniform" kernel. Other options are "Epanichnekov" and "normal".
#'        It will only be used when the bandwidth \code{h} is specified.
#' @param x0 a scalar as the center around which the kernel weights are calculated.
#' @param X the vector of biomarker values from other subjects used to calculate the weights around the center x0.
#' @return Return a vector of kernel weights for each element in X. It has the same length as X.
#' @note X must be the vector of ALL biomarker values in the data; it cannot
#'       be any other vector of arbitrary length.
#' @keywords internal


calc.kw <- function( X, x0, span = 0.1, h=NULL, type="uniform" ) {
  # Calculate the nearest neighbor kernel weights
  # Arguments:
  #  -- span: the proportion of observations used
  #  -- h: the bandwidth of kernel weights
  #  -- type: the type of kernel function, "uniform", "Epanichnekov" and "normal"
  #  -- x0: the x0 around which the kernel weights are calculated
  #  -- X: the vector of all biomarker values in the data
  # Return:
  #  a vector of kernel weights for each element in X
  # NOTE: X must be the vector of ALL biomarker values in the data; it cannot
  #       be any other vector of arbitrary length

  if (is.null (h) & !is.null(span)){
    n <- length(X) ;
    tmp0 <- abs( X-x0 ) ;
    tmp1 <- sort( tmp0 ) ;
    tmp2 <- tmp1[ ceiling(n*span) ] ;
    # the cut off that defines the neighborhood
    ans <- as.numeric( tmp0 <= tmp2 ) ;
  } else if (!is.null(h)){
    x.new <- (X-x0)/h ;

    if ( type == "uniform" ) {
      ans <- as.numeric( abs(x.new) <= 1 )/(2*h) ;
    } else if ( type == "Epanichnekov" ) {
      ans <- 0.75*(1-x.new*x.new)*as.numeric( abs(x.new) <= 1 )/h ;
    } else {
      ans <- 1/sqrt(2*pi)*exp(-1/2*x.new^2)/h
    }
  } else if ( is.null(span) & is.null(h)){
    print ("error! span or h needs to be specified!")
  }

  ans ;
}

