#' The full Fisher Information Matrix (FIM) Calculating one model switch at a time, good for large matrices.
#' 
#' Compute the full FIM given specific model(s), parameters, design and methods. 
#' This computation calculates the FIM for each model switch separately.  Correlations between the models parameters are assumed to be zero.
#' 
#' @inheritParams mftot
#' 
#' @return As a list:
#' \item{ret}{The FIM}
#' \item{poped.db}{A PopED database}
#' 
#' @seealso For an easier function to use, please see \code{\link{evaluate.fim}}.  
#' @family FIM
#' 
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_mftot6.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mftot6 <- function(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db){
m=size(ni,1)
s=0
for(i in 1:m){
   if((ni[i]!=0 && groupsize[i]!=0)){
      
         if((!isempty(x))){
           x_i = t(x[i,,drop=F])      
         } else {
           x_i =  zeros(0,1)
         }
         if((!isempty(a))){
           a_i = t(a[i,,drop=F])
         } else {
           a_i =  zeros(0,1)
         }
      
         returnArgs <- mf7(t(model_switch[i,1:ni[i,drop=F],drop=F]),t(xt[i,1:ni[i,drop=F],drop=F]),x_i,a_i,bpop,d,sigma,docc,poped.db) 
mf_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
        s=s+groupsize[i]*mf_tmp
    }
}
ret = s
return(list( ret= ret,poped.db =poped.db )) 
}
