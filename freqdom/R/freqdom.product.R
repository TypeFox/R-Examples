#' For given spectral densities \code{S} and \code{SC} computes product \code{S \%*\% SC} at each frequency.
#' 
#' @title Compute a product of two spectral densities 
#' @param S first spectral density
#' @param SC second spectral density
#' @return Frequency Domain Operator object
#' @seealso \code{\link{freqdom.inverse}}, \code{\link{freqdom.ratio}}
#' @export
#' @examples
#' n = 100
#' X = rar(n)
#' Y = rar(n)
#' SX = spectral.density(X)
#' SY = spectral.density(Y)
#' R = freqdom.product(SY,SX)
freqdom.product = function(S,SC){
  if (!is.freqdom(S))
    stop("S must be a freqdom object")
  if (!is.freqdom(SC))
    stop("SC must be a freqdom object")
  if (dim(S$operators)[3] != dim(SC$operators)[2])
    stop("Dimensions of operators don't match")
  
  R = SC
  D = c(dim(S$operators)[1],
        dim(S$operators)[2],
        dim(SC$operators)[3])
  R$operators = array(0,D)
  for (theta in 1:length(S$freq))
    R$operators[theta,,] = S$operators[theta,,] %*% SC$operators[theta,,]
  R
}

oldprod <- `%*%`

#' Frequency-wise or time-wise matrix product. Takes two elements
#' \code{freqdom} or \code{timedom} and multiplies them on
#' each frequency or time point. If objects of other type are provided
#' then the standard multiplication is applied.
#'  
#' @title Frequency-wise or component-wise matrix product. 
#' @param e1 first element
#' @param e2 second element
#' @return object of the same type as e1 but with new dimensions
#' @export
`%*%` <- function (e1,e2) {
  if (is.freqdom(e1) && is.freqdom(e2))
    freqdom.product(e1,e2)
  else
    oldprod(e1,e2)
}
