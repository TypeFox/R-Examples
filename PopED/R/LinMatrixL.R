#' The linearized matrix L
#' 
#' Function computes the derivative of the model with respect to the between subject variability 
#' terms in the model (b's and bocc's) evaluated at
#' a defined point 
#' (b_ind and bocc_ind).
#' 
#' @inheritParams mf
#' @param bpop The fixed effects parameter values.  Supplied as a vector.
#' @param b_ind The point at which to evaluate the derivative
#' @param bocc_ind The point at which to evaluate the derivative
#' @param poped.db A PopED database.
#' 
#' @return As a list:
#' \item{y}{A matrix of size (samples per individual x number of random effects)}
#' \item{poped.db}{A PopED database}
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_LinMatrixL.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

LinMatrixL <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db){
  
  if((poped.db$parameters$NumRanEff==0)){
    y=0
  } else {
    returnArgs <- gradff(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db) 
    grad_ff_tmp <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    y=grad_ff_tmp%*%gradfg(x,a,bpop,b_ind,bocc_ind,poped.db)
  }
  return(list( y= y,poped.db=poped.db)) 
}


