#' The reduced Fisher Information Matrix (FIM) parameterized with A,B,C matrices & using the derivative of variance.
#' 
#' Compute the reduced FIM given specific model(s), parameters, design and methods. 
#' This computation assumes that there is no correlation in the FIM between the fixed and random effects, 
#' and set these elements in the FIM to zero.
#' This computation parameterizes the FIM calculation using 
#' A,B,C matrices (as in Retout \emph{et al.}) but uses the derivative of variances.
#' Should give the same answer as \code{\link{mftot1}} but computation times may be different.   
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
#' @references S. Retout and F. Mentre, "Further developments of the Fisher Information Matrix in
#' nonlinear mixed effects models with evaluation in population pharmacokinetics", J. of Biopharm. Stats., 13(2), 2003.
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_mftot7.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

mftot7 <- function(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db){
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

         returnArgs <- mf8(t(model_switch[i,1:ni[i,drop=F],drop=F]),t(xt[i,1:ni[i,drop=F],drop=F]),x_i,a_i,bpop,d,sigma,docc,poped.db) 
mf_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
        s=s+groupsize[i]*mf_tmp
    }
}
ret = s
return(list( ret= ret,poped.db =poped.db )) 
}
