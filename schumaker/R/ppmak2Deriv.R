#' ppmak2Deriv
#'
#' Create the second derivative of the spline defined by given intervals and quadratic coefficients.
#' This is an internal function that is called from the Schumaker function.
#' @param IntStarts This is a vector with the start of each interval.
#' @param SpCoefs This is a matrix with three columns. The first is the coefficient of the squared term followed by linear term coefficients and constants.
#' @param Vectorised This is a boolean parameter. Set to TRUE if you want to be able to input vectors to the created spline. If you will only input single values set this to FALSE as it is a bit faster.

#' @return A spline function for the given intervals and quadratic curves. Each function takes an x value (or vector if Vectorised = TRUE) and outputs the interpolated y value (or relevent derivative).
ppmak2Deriv = function(IntStarts, SpCoefs, Vectorised  = TRUE){
  if (!(Vectorised)){
    sp = function(PointToExamine){
      IntervalNum = findInterval(PointToExamine, IntStarts, all.inside = TRUE)
      xmt = PointToExamine - IntStarts[IntervalNum]
      Coefs = SpCoefs[IntervalNum,,drop = FALSE]
      Coefs %*% c(2,0, 0)
    }
  } else {
    sp = function(PointToExamine){
      IntervalNum = findInterval(PointToExamine, IntStarts, all.inside = TRUE)
      xmt = PointToExamine - IntStarts[IntervalNum]
      Len = length(PointToExamine)
      xmtMat = matrix(c(rep(2, Len),rep(0, Len), rep(0, Len)), ncol = 3, byrow = FALSE)
      Coefs = SpCoefs[IntervalNum,,drop = FALSE]
      (Coefs * xmtMat) %*% c(1,1,1)
    }
  }
}
