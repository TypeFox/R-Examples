#' Downsize a general design to a specific design
#' 
#' Function takes a design with potentially empty design 
#' variables and resuces the design so that a FIM can be calculated using \code{\link{mftot}}.
#' 
#' @param poped.db A PopED database 
#' @return A list containing:
#' \item{ni}{A vector of the number of samples in each group.}
#' \item{xt}{A matrix of sample times.  Each row is a vector of sample times for a group.}
#' \item{model_switch}{A matrix that is the same size as xt, specifying which model each sample belongs to.}
#' \item{x}{A matrix for the discrete design variables.  Each row is a group.}
#' \item{a}{A matrix of covariates.  Each row is a group.}
#' \item{bpop}{A matrix of fixed effect parameter values.}
#' 
#' @family poped_input
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_downsizing_general_design.R
#' @export
#' @keywords internal
#' 
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker



downsizing_general_design <- function(poped.db){
  # ------------- downsizing of general design
  
  ni=poped.db$design$ni[1:poped.db$design$m,,drop=F]
  xt=poped.db$design$xt[1:poped.db$design$m,1:max(poped.db$design_space$maxni),drop=F]
  model_switch = poped.db$design$model_switch[1:poped.db$design$m,1:max(poped.db$design_space$maxni),drop=F]
  
  if((size(poped.db$design$x,2)!=0)){
    x=poped.db$design$x[1:poped.db$design$m,1:size(poped.db$design$x,2),drop=F]
  } else {
    x = zeros(poped.db$design$m,0)
  }
  if((size(poped.db$design$a,2)!=0)){
    a=poped.db$design$a[1:poped.db$design$m,1:size(poped.db$design$a,2),drop=F]
    maxa=poped.db$design_space$maxa[1:poped.db$design$m,1:size(poped.db$design$a,2),drop=F]
    mina=poped.db$design_space$mina[1:poped.db$design$m,1:size(poped.db$design$a,2),drop=F]
  } else {
    a = zeros(poped.db$design$m,0)
    maxa = matrix(0,0,0)
    mina = matrix(0,0,0)
  }
  bpop=poped.db$parameters$bpop[1:poped.db$parameters$nbpop,1:3,drop=F]
  n=t(ni)%*%matrix(1,poped.db$design$m,1)
  
  d=poped.db$parameters$d[1:poped.db$parameters$NumRanEff,1:3,drop=F]
  maxxt=poped.db$design_space$maxxt[1:poped.db$design$m,1:max(poped.db$design_space$maxni),drop=F]
  minxt=poped.db$design_space$minxt[1:poped.db$design$m,1:max(poped.db$design_space$maxni),drop=F]
  
  return(list( ni= ni, xt= xt, model_switch= model_switch, x= x, a= a, bpop= bpop, 
               n= n, d= d, maxxt= maxxt, minxt= minxt,maxa=maxa,mina =mina )) 
}
