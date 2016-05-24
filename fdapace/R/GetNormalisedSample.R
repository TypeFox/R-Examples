#' Normalise sparse functional sample
#'
#' Normalise sparse functional sample given in an FPCA object
#'
#' @param fpcaObj An FPCA object.
#'
#' @return A list containing the normalised sample 'y' at times 't'
#'
#' @export
GetNormalisedSample<- function(fpcaObj){
  if (any( 0>=diag(fpcaObj$smoothedCov)) ){
    stop("The autocovariance functions appears to have negative or zero values.")
  }
  ynorm = mapply(FUN = function(vy, vt){ 
    return( ( vy - approx(y = fpcaObj$mu, x =fpcaObj$workGrid, xout = vt)$y) /
              approx(y = sqrt(diag(fpcaObj$smoothedCov)), x =fpcaObj$workGrid, xout = vt)$y)
  }, vy = fpcaObj$inputData$y, vt = fpcaObj$inputData$t, SIMPLIFY = FALSE)
  return(list(y = ynorm, t = fpcaObj$inputData$t ))
}
