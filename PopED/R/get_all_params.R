#' Extract all model parameters from the PopED database.
#' 
#' @param poped.db A PopED database.
#' @return A list containing:
#' \item{bpop}{A vector of fixed effect parameter values.}
#' \item{d}{A vector of between subject variability parameters}
#' \item{covd}{A vector of the covariances of the between subject variability parameters.  Row major format of the lower triangular portion of the D (OMEGA) matrix}
#' \item{docc}{A vector of the between occasion variability (BOV) terms in the model}
#' \item{covdocc}{A vector of the covariances between the BOV terms.  Row major of the lower triangular portion of the BOV matrix. }
#' \item{sigma}{A vector of the resudual unexplained variances (RUV)}
#' \item{covsigma}{A vector of the covariances between the RUV terms}
#' \item{all}{A vector with all of the above, in the order of this list.}
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_get_all_params.R
#' @export
#' @keywords internal

## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

get_all_params <- function(poped.db){
  #Return all params (in a vector all) with the specified order above
  
  bpop = poped.db$parameters$bpop[,2,drop=F]
  d = poped.db$parameters$d[,2,drop=F]
  docc = poped.db$parameters$docc[,2,drop=F]
  covd = poped.db$parameters$covd
  sigma = diag_matlab(poped.db$parameters$sigma)
  covsigma = zeros(1,length(sigma)*(length(sigma)-1)/2)
  
  k=1
  for(i in 1:size(poped.db$parameters$sigma,1)){
    for(j in 1:size(poped.db$parameters$sigma,2)){
      if((i<j)){
        covsigma[k] = poped.db$parameters$sigma[i,j]
        k=k+1
      }
    }
  }
  
  covdocc = poped.db$parameters$covdocc
  
  all = matrix(c(bpop, d, t(covd), docc, t(covdocc), sigma, t(covsigma)),ncol=1,byrow=T)
  return(list( bpop= bpop,d=d,covd=covd,docc=docc,covdocc=covdocc,sigma=sigma,covsigma=covsigma,all =all )) 
}
