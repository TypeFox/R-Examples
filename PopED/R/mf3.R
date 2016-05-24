#' The reduced  Fisher Information Matrix (FIM) for one individual
#' 
#' Compute the reduced  FIM for one individual given specific model(s), parameters, design and methods. 
#' This computation assumes that there is no correlation in the FIM between the fixed and random effects, 
#' and set these elements in the FIM to zero.
#' 
#' @param xt A vector of sample times.  
#' @inheritParams mf
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{poped.db}{A PopED database}
#' 
#' @seealso Used by \code{\link{mftot1}}.  
#' @family FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_mf3.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf3 <- function(model_switch,xt,x,a,bpop,d,sigma,docc,poped.db){
  #Calculate the reduced FIM
  
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
    
    if((poped.db$settings$bCalculateEBE) ){#Calculate an EBE
      epsi0 = zeros(1,length(poped.db$parameters$notfixed_sigma))
      g=feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind)
      returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt,g,epsi0,poped.db) 
      mean_data <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      start_bind = t(b_ind)
      b_ind = ind_estimates(mean_data,bpop,d,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt,x,a,b_ind,bocc_ind,poped.db)
      #        b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(b_ind),(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
      
      #b_ind2 = ind_estimates(mean_data,bpop,d,sigma,t(zeros(size(b_ind)[1],size(b_ind)[2])),!(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
      poped.db$mean_data = mean_data
    }
    
    n_fixed_eff <- numnotfixed_bpop
    n_rand_eff <- numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma
    start_col <- 1
    
    f1=zeros(n+n*n,n_fixed_eff+n_rand_eff)
    
    if(n_fixed_eff!=0){
      returnArgs <- m1(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,poped.db) 
      f1[1:n,start_col:(start_col+n_fixed_eff-1)] <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      start_col <- start_col + n_fixed_eff
    }
    if(n_rand_eff!=0){
      returnArgs <- m3(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,TRUE,poped.db) 
      f1[(n+1):(n+n*n),start_col:(start_col+n_rand_eff-1)] <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
    }
    
    f2=zeros(n+n*n,n+n*n)
    returnArgs <-  v(model_switch,xt,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
    v_tmp <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    if((matrix_any(v_tmp)!=0) ){#If the inverse is not empty
      f2[1:n,1:n]=inv(v_tmp)
      tmp_m4=m4(v_tmp,n)
      f2[(n+1):(n+n*n),(n+1):(n+n*n)]=inv(tmp_m4)
    }
    if((all(f2==0))){
      ret = ret+t(f1)*f1
    } else {
      ret = ret+t(f1)%*%f2%*%f1
    }
  }
  ret = ret/poped.db$settings$iFOCENumInd
  
  return(list( ret= ret,poped.db=poped.db)) 
}
