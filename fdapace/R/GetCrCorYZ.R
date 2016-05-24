#' Make cross-correlation matrix from auto- and cross-covariance matrix
#'
#' Make cross-correlation matrix from auto- andcross-covariance matrix
#'
#' @param ccYZ The cross-covariance vector between variables Y and Z (n-by-1).
#' @param acYY The auto-covariance n-by-n matrix of variable Y or the (n-by-1) diagonal of that matrix.
#' @param covZ The (scalar) covariance of variable Z. 
#' 
#' @return A cross-correlation matrix between variables Y (functional) and Z (scalar).
#'
#' @export
#' 
GetCrCorYZ <- function(ccYZ, acYY , covZ){
  
  acYY = as.matrix(acYY); # Such a messy things because R does not treat vectors as matrices.
  
  # Check basic sizes
  if(1 != length(covZ)){
    stop('The variance of Z must be a scalar.')
  }
  if(!is.matrix(ccYZ) && !is.vector(ccYZ)){
    stop('The cross-covariance of Y and Z is must be a matrix or a vector.')
  }
  if( is.matrix(ccYZ) && (1 != dim(ccYZ)[2]) ){
    stop('The cross-covariance matrix of a functional variable Y and a scalar variable Z is not n-by-1.')
  }
  N = length(ccYZ);
  
  if( N^2 == length(acYY)){
    diagYY = diag(acYY)
  } else {
    if( 1 != dim(acYY)[2]){
      stop('The auto-covariance is not n-by-n or n-by-1 ')
    } else {
      diagYY = as.vector(acYY)
    }
  }
  
  diagZ =(covZ[1])
  
  if( length(diagYY) != length(ccYZ) ){
    stop('The cross-covariance for YZ and the provided covariance for Y are of incompatible sizes.')
  }
  if( any(1e-12> diagYY)){
    stop('The provided covariance for X are unreasonable small or negative. Rescale/check your data.')
  }
  if( any(1e-12> diagZ)){
    stop('The provided covariance for Z are unreasonable small or negative. Rescale/check your data.')
  }
  
  # return (solve(sqrt(diag(diagYY))) %*% ccYZ %*% solve(sqrt(diag(diagZ))))
  return( diag(1/sqrt(diagYY) , nrow = length(diagYY)) %*% ccYZ %*% diag(1/sqrt(diagZ), nrow = length(diagZ)) )
}



