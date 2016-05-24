#' Create a full D (between subject variability) matrix given a vector of variances and covariances.
#' 
#' @param variance_vector The vector of the variances.
#' @param covariance_vector A vector of the covariances. Writen in column major 
#'   order for the lower triangular matrix.
#' @return The full matrix of variances for the between subject variances
#' @example tests/testthat/examples_fcn_doc/examples_getfulld.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

getfulld <- function(variance_vector,covariance_vector=NULL){  
  if(length(variance_vector)==1) return(as.matrix(variance_vector))
  if((isempty(covariance_vector) || sum(covariance_vector!=0)==0)){
    d=diag_matlab(variance_vector)
  } else {
    covd = zeros(length(variance_vector),length(variance_vector))
    k=1
    for(i in 1:length(variance_vector)){
      for(j in 1:length(variance_vector)){
        if((i<j)){
          covd[i,j]=covariance_vector[k]
          covd[j,i]=covariance_vector[k]
          k=k+1
        }
      }
    }
    d = diag_matlab(variance_vector)+covd
  }
  return(d) 
}
