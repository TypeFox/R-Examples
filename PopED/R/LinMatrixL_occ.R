#' Model linearization with respect to occasion variablity parameters.
#' 
#' The function performs a linearization of the model with respect to the occation  variability parameter..
#' Derivative of model w.r.t. eta_occ, evaluated bocc_ind.
#' 
#' @inheritParams mftot
#' @inheritParams LinMatrixH
#' @param iCurrentOcc The current occasion.
#' 
#' @return A matrix of size (samples per individual x number of iovs)
#'  
#' @family FIM
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_LinMatrixL_occ.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

LinMatrixL_occ <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,iCurrentOcc,poped.db){
#
# size: (samples per individual x number of iovs)
#
if((poped.db$parameters$NumOcc==0)){
	y=0
} else {
     returnArgs <- gradff(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db) 
grad_ff_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
    y=grad_ff_tmp%*%gradfg_occ(x,a,bpop,b_ind,bocc_ind,iCurrentOcc,poped.db)
}
return(list( y= y,poped.db=poped.db)) 
}


