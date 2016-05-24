#' The reduced Fisher Information Matrix (FIM) for one individual parameterized with A,B,C matrices & using the derivative of variance.
#' 
#' Compute the reduced FIM for one individual given specific model(s), parameters, design and methods. 
#' This computation assumes that there is no correlation in the FIM between the fixed and random effects, 
#' and set these elements in the FIM to zero.
#' This computation parameterizes the FIM calculation using 
#' A,B,C matrices (as in Retout \emph{et al.}) but uses the derivative of variances.
#' Should give the same answer as \code{\link{mf3}} but computation times may be different.   
#' 
#' @inheritParams mf
#' 
#' @return As a list:
#' \item{ret}{The FIM for one individual}
#' \item{poped.db}{A PopED database}
#' 
#' @seealso Used by \code{\link{mftot7}}.  
#' @family FIM
#' 
#' @references S. Retout and F. Mentre, "Further developments of the Fisher Information Matrix in
#' nonlinear mixed effects models with evaluation in population pharmacokinetics", J. of Biopharm. Stats., 13(2), 2003.
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_mf8.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mf8 <- function(model_switch,xt_ind,x,a,bpop,d,sigma,docc,poped.db){

#Calculate FIM with another parameterization, i$e. the parametrization used
#in Retout et al but with derivative of variance instead of SD for sigma and the reduced FIM

numnotfixed_bpop = sum(poped.db$parameters$notfixed_bpop)
numnotfixed_d    = sum(poped.db$parameters$notfixed_d)
numnotfixed_covd = sum(poped.db$parameters$notfixed_covd)
numnotfixed_docc  = sum(poped.db$parameters$notfixed_docc)
numnotfixed_covdocc  = sum(poped.db$parameters$notfixed_covdocc)
numnotfixed_sigma  = sum(poped.db$parameters$notfixed_sigma)
numnotfixed_covsigma  = sum(poped.db$parameters$notfixed_covsigma)

n=size(xt_ind,1)
ret = 0

for(i in 1:poped.db$settings$iFOCENumInd){
    b_ind = poped.db$parameters$b_global[,i,drop=F]
    bocc_ind = poped.db$parameters$bocc_global[[i]]
    tmp_fim=zeros(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma,numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)
     returnArgs <-  m1(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,poped.db) 
m1_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
     returnArgs <-  m3(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,TRUE,poped.db) 
m3_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
     returnArgs <-  v(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
v_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
    invv = inv(v_tmp)
    
    tmp_fim[1:numnotfixed_bpop,1:numnotfixed_bpop]=2*t(m1_tmp)%*%invv%*%m1_tmp
    for(m in 1:numnotfixed_bpop){
        for(k in 1:numnotfixed_bpop){
            tmp_fim[m,k]=tmp_fim[m,k]
        }
    }

    for(m in 1:(numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)){
        for(k in 1:(numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)){
            tmp_fim[numnotfixed_bpop+m,numnotfixed_bpop+k]=trace_matrix(reshape_matlab(m3_tmp[,m,drop=F],n,n)%*%invv%*%reshape_matlab(m3_tmp[,k,drop=F],n,n)%*%invv)
        }
    }

    ret = ret+1/2*tmp_fim
}
ret = ret/poped.db$settings$iFOCENumInd
return(list( ret= ret,poped.db=poped.db)) 
}

