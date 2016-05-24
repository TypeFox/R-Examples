#' For given spectral densities \code{S} and \code{SC} compute the Kronecker product
#' \code{S \%x\% SC} at each frequency.
#' 
#' @title Compute a kronecker product of two spectral densities 
#' @param S first spectral density
#' @param SC second spectral density
#' @return Frequency Domain Operator object
#' @export
#' @examples
#' n = 100
#' X = rar(n,d=3)
#' Y = rar(n,d=3)
#' SX = spectral.density(X)
#' SY = spectral.density(Y)
#' R = freqdom.kronecker(SY,SX)
freqdom.kronecker = function(S,SC){
  if (!is.freqdom(S))
    stop("S must be a freqdom object")
  if (!is.freqdom(SC))
    stop("SC must be a freqdom object")
  
  R = SC

  D1 = dim(S$operators)
  D2 = dim(SC$operators)
  D = c(0,0,0)
  D[2] = D1[2]*D2[2]
  D[3] = D1[3]*D2[3]
  D[1] = D1[1]
  R$operators = array(0,D)

  for (theta in 1:length(S$freq))
    R$operators[theta,,] = S$operators[theta,,] %x% SC$operators[theta,,]
  R
}

oldkronprod <- `%x%`

#' Frequency-wise or time-wise Kronecker product. Takes two elements
#' \code{freqdom} or \code{timedom} and applies the Kronecker product on
#' each frequency or time point. If objects of other type are provided
#' then the standard function is applied.
#'  
#' @title Frequency-wise or component-wise Kronecker product. 
#' @param e1 first element
#' @param e2 second element
#' @return object of the same type as e1 but with new dimensions
#' @export
`%x%` <- function (e1,e2) {
  if (is.freqdom(e1) && is.freqdom(e2))
    freqdom.kronecker(e1,e2)
  else
    oldkronprod(e1,e2)
}
