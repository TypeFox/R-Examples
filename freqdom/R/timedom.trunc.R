#' Truncates a time domain object to specified lags (equivalent to setting the lags to zero)
#'
#' @title Truncates a time domain object to specified lags
#' @param A \code{timedom} object - series of operators
#' @param lags lags to which A should be truncated
#' @return Truncated time series
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' #estimate regressors in model $Y_t = \sum_{i\in Z} A_i X_{t-i}$
#' A = speclagreg(X,Y)
#' B = timedom.trunc(A, c(-1, 2, 3))
#' @export 
timedom.trunc = function(A, lags){
  if (!is.timedom(A))
    stop ("A must be a time domain filter")
  
  D = dim(A$operators)
  D[1] = sum(A$lags %in% lags)
  A$operators = array(A$operators[A$lags %in% lags,,],D)
  A$lags = intersect(lags,A$lags)
  A
}

