#' The full Fisher Information Matrix (FIM) for one individual Calculating one model switch at a time, good for large matrices.
#' 
#' Compute the full FIM for one individual given specific model(s), parameters, design and methods. 
#' This computation calculates the FIM for each model switch separately.  Correlations between the models parameters are assumed to be zero.
#' 
#' @inheritParams mf
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{poped.db}{A PopED database}
#' 
#' @seealso Used by \code{\link{mftot6}}.  
#' @family FIM
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_mf7.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf7 <- function(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped.db){

#This calculation of FIM devide the calculation up into one calculation
#per model switch

numnotfixed_bpop = sum(poped.db$parameters$notfixed_bpop)
numnotfixed_d    = sum(poped.db$parameters$notfixed_d)
numnotfixed_covd = sum(poped.db$parameters$notfixed_covd)
numnotfixed_docc  = sum(poped.db$parameters$notfixed_docc)
numnotfixed_covdocc  = sum(poped.db$parameters$notfixed_covdocc)
numnotfixed_sigma  = sum(poped.db$parameters$notfixed_sigma)
numnotfixed_covsigma  = sum(poped.db$parameters$notfixed_covsigma)

ret = 0

for(i in 1:poped.db$settings$iFOCENumInd){
    b_ind = poped.db$parameters$b_global[,i,drop=F]
    bocc_ind = poped.db$parameters$bocc_global[[i]]
    
    for(j in 1:max(model_switch)){
        xt_new = zeros(sum(model_switch==j),1)
        m=1
        for(k in 1:length(model_switch)){
            if((model_switch[k]==j)){
                xt_new[m]=xt_ind[k]
                m=m+1
            }
        }
        n=size(xt_new,1)
        model_switch_new = matrix(1,sum(model_switch==j),1)*j
        
        if((!isempty(xt_new))){
            f1=zeros(n+n*n,numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)
             returnArgs <- m1(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,poped.db) 
f1[1:n,1:numnotfixed_bpop] <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
             returnArgs <- m2(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
f1[(n+1):(n+n*n),1:numnotfixed_bpop] <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
             returnArgs <- m3(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,TRUE,poped.db) 
f1[(n+1):(n+n*n),(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)] <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
            f2=zeros(n+n*n,n+n*n)
             returnArgs <-  v(model_switch_new,xt_new,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
v_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
            if((matrix_any(v_tmp)!=0)){
               # f2[1:n,1:n]=v_tmp\diag_matlab(1,size(v_tmp))
                f2[1:n,1:n]=inv(v_tmp)
                #m4_tmp = m4(v_tmp,n)
                #f2[(n+1):(n+n*n),(n+1):(n+n*n)]=m4_tmp\diag_matlab(1,size(m4_tmp))
                f2[(n+1):(n+n*n),(n+1):(n+n*n)]=inv(m4(v_tmp,n))
            }
            ret=ret+t(f1)%*%f2%*%f1
        }
    }
}
ret = ret/poped.db$settings$iFOCENumInd
return(list( ret= ret,poped.db=poped.db)) 
}

