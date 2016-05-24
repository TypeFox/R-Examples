#' The reduced Fisher Information Matrix (FIM) for one individual, using the SD of RUV as a parameter.
#'  
#' 
#' Compute the reduced FIM for one individual using the standard deviation of the residual unexplained variability (RUV) terms as a parameter, 
#' given specific model(s), parameters, design and methods. 
#' This computation 
#' assumes that there is no correlation in the FIM between the fixed and random effects, 
#' and set these elements in the FIM to zero.
#' In addition all derivatives in the computation are made 
#' with respect to the standard deviation of the RUV terms (sqrt(SIGMA) in NONMEM). 
#' This matches what is done in PFIM, and assumes that the standard deviation of the residual unexplained variation is the estimated parameter
#' (NOTE: NONMEM estimates the variance of the resudual unexplained variation by default).
#' 
#' @inheritParams mf3
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{poped.db}{A PopED database}
#' 
#' @seealso Used by \code{\link{mftot4}}.  
#' @family FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_mf5.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf5 <- function(model_switch,xt,x,a,bpop,d,sigma,docc,poped.db){
  #==== Reduced FIM with derivative of variance w$r.t. the SD, NOT the
  # variance
  
  
  numnotfixed_bpop = sum(poped.db$parameters$notfixed_bpop)
  numnotfixed_d    = sum(poped.db$parameters$notfixed_d)
  numnotfixed_covd = sum(poped.db$parameters$notfixed_covd)
  numnotfixed_docc  = sum(poped.db$parameters$notfixed_docc)
  numnotfixed_covdocc  = sum(poped.db$parameters$notfixed_covdocc)
  numnotfixed_sigma  = sum(poped.db$parameters$notfixed_sigma)
  numnotfixed_covsigma  = sum(poped.db$parameters$notfixed_covsigma)
  
  
  n=size(xt,1)
  ret = 0
  
  for(i in 1:poped.db$settings$iFOCENumInd){
    b_ind = poped.db$parameters$b_global[,i,drop=F]
    bocc_ind = poped.db$parameters$bocc_global[[i]]
    f1=zeros(n+n*n,numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)
    returnArgs <- m1(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,poped.db) 
    f1[1:n,1:numnotfixed_bpop] <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    #w$r.t. the variance
    #[f1[(n+1):(n+n*n),(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)],poped.db]=m3(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,TRUE,poped.db)
    #w$r.t. the sd
    returnArgs <- m3(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,FALSE,poped.db) 
    f1[(n+1):(n+n*n),(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)] <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    
    f2=zeros(n+n*n,n+n*n)
    returnArgs <-  v(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
    v_tmp <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    if((matrix_any(v_tmp)!=0) ){#If the inverse is not empty
      f2[1:n,1:n]=inv(v_tmp)
      tmp_m4=m4(v_tmp,n)
      f2[(n+1):(n+n*n),(n+1):(n+n*n)]=inv(tmp_m4)
    }
    ret = ret+t(f1)%*%f2%*%f1
  }
  ret = ret/poped.db$settings$iFOCENumInd
  return(list( ret= ret,poped.db=poped.db)) 
}

