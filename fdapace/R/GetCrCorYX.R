#' Make cross-correlation matrix from auto- and cross-covariance matrix
#'
#' Make cross-correlation matrix from  auto- andcross-covariance matrix
#'
#' @param ccXY The cross-covariance matrix between variables X and Y.
#' @param ccXX The auto-covariance matrix of variable X or the diagonal of that matrix.
#' @param ccYY The auto-covariance matrix of variable Y or the diagonal of that matrix.
#' 
#' @return A cross-correlation matrix between variables X and Y.
#'
#' @export
#' 
#' 
GetCrCorYX <- function(ccXY, ccXX , ccYY){
  
  if(!is.matrix(ccXY)){
    stop('The cross-covariance matrix is must be a matrix.')
  }
  
  if(!is.matrix(ccXX) && !is.vector(ccXX)){
    stop('The auto-covariance matrix for X must be a matrix or vector.')
  }
  
  if(!is.matrix(ccYY) && !is.vector(ccYY)){
    stop('The auto-covariance matrix for Y must be a matrix or vector.')
  }
  
  if(is.matrix(ccYY)){
    diagYY = diag(ccYY)
  } else {
    diagYY =(ccYY)
  }  
  
  if(is.matrix(ccXX)){
    diagXX = diag(ccXX)
  } else {
    diagXX =(ccXX)
  }  
  
  if( length(diagXX) != dim(ccXY)[1] ){
    stop('The cross-covariance matrix for XY and the provided covariance for X are incompatible.')
  }
  
  if( length(diagYY) != dim(ccXY)[2] ){
    stop('The cross-covariance matrix for XY and the provided covariance for Y are incompatible.')
  }
  
  if( any(1e-12> diagXX)){
    stop('The provided covariance for X are unreasonable small or negative. Rescale/check your data.')
  }
  
  if( any(1e-12> diagYY)){
    stop('The provided covariance for X are unreasonable small or negative. Rescale/check your data.')
  }
  
  # return (solve(sqrt(diag(diagXX))) %*% ccXY %*% solve(sqrt(diag(diagYY))))
  return( diag(1/sqrt(diagXX) , nrow = length(diagXX)) %*% ccXY %*% diag(1/sqrt(diagYY), nrow = length(diagYY)) )
}


 
