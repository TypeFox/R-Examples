#' Symmetric square root of a matrix
#' 
#' Computes the symmetric square root of a positive definite matrix
#' 
#' 
#' @usage mhalf(M)
#' @param M a positive definite matrix
#' @return a matrix \code{H} such that \code{H^2} equals \code{M}
#' @author Peter Hoff
#' @export mhalf
mhalf <-
function(M) 
{ 
  #symmetric square  root of a pos def matrix
  tmp<-eigen(M)
  tmp$vec%*%sqrt(diag(tmp$val,nrow=nrow(M)))%*%t(tmp$vec)
}
