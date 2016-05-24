#' Evaluate the Fisher Information Matrix (FIM)
#' 
#' Compute the FIM given the model, parameters, design and methods defined in the 
#' PopED database. Some of the arguments coming from the PopED database can be overwritten;  
#' by default these arguments are \code{NULL} in the 
#' function, if they are supplied then they are used instead of the arguments from the PopED database.
#' 
#' @param poped.db A PopED database.
#' @param fim.calc.type The method used for calculating the FIM. Potential values:
#' \itemize{
#' \item 0 = Full FIM.  No assumption that fixed and random effects are uncorrelated.  See \code{\link{mftot0}}.
#' \item 1 = Reduced FIM. Assume that there is no correlation in the FIM between the fixed and random effects, and set these elements in 
#' the FIM to zero. See \code{\link{mftot1}}.
#' \item 2 = weighted models (placeholder).
#' \item 3 = Not currently used.
#' \item 4 = Reduced FIM and computing all derivatives with respect to the standard deviation of the residual unexplained variation (sqrt(SIGMA) in NONMEM). 
#' This matches what is done in PFIM, and assumes that the standard deviation of the residual unexplained variation is the estimated parameter
#' (NOTE: NONMEM estimates the variance of the resudual unexplained variation by default). See \code{\link{mftot4}}.
#' \item 5 = Full FIM parameterized with A,B,C matrices & derivative of variance. See \code{\link{mftot5}}.
#' \item 6 = Calculate one model switch at a time, good for large matrices. See \code{\link{mftot6}}.
#' \item 7 = Reduced FIM parameterized with A,B,C matrices & derivative of variance See \code{\link{mftot7}}.
#' }
#' @param approx.method Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI
#' @param FOCE.num Number indivduals in each step of FOCE approximation method 
#' @param bpop.val The fixed effects parameter values.  Supplied as a vector.
#' @param d_full A between subject variability matrix (OMEGA in NONMEM).
#' @param docc_full A between occasion variability matrix.
#' @param sigma_full A residual unexplained variability matrix (SIGMA in NONMEM).
#' @param model_switch A matrix that is the same size as xt, specifying which model each sample belongs to.
#' @param ni A vector of the number of samples in each group.
#' @param xt A matrix of sample times.  Each row is a vector of sample times for a group.
#' @param x A matrix for the discrete design variables.  Each row is a group.
#' @param a A matrix of covariates.  Each row is a group.
#' @param groupsize A vector of the numer of individuals in each group.
#' @param deriv.type A number indicating the type of derivative to use:
#' \itemize{
#' \item 0=Complex difference 
#' \item 1=Central difference 
#' \item 20=Analytic derivative (placeholder) 
#' \item 30=Automatic differentiation (placeholder)
#' }
#' @param ... Other arguments passed to the function.
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' 
#' @return The FIM.
#' 
#' @family FIM
#' @family evaluate_design
#' @family evaluate_FIM
#' 
# @example inst/examples_fcn_doc/examples_evaluate.fim.R
#' 
#' @example tests/testthat/examples_fcn_doc/examples_evaluate.fim.R
#' @export


evaluate.fim <- function(poped.db,
                         fim.calc.type=NULL,
                         approx.method=NULL, 
                         FOCE.num = NULL,
                         bpop.val=NULL,
                         d_full=NULL,
                         docc_full=NULL,
                         sigma_full=NULL,
                         model_switch=NULL,
                         ni=NULL,
                         xt=NULL,
                         x=NULL,
                         a=NULL,
                         groupsize=NULL,
                         deriv.type = NULL,
                         ...){
  
  
  if(is.null(bpop.val)) bpop.val <- poped.db$parameters$param.pt.val$bpop
  if(is.null(d_full)) d_full <- poped.db$parameters$param.pt.val$d
  if(is.null(docc_full)) docc_full <- poped.db$parameters$param.pt.val$docc
  if(is.null(sigma_full)) sigma_full <- poped.db$parameters$param.pt.val$sigma
  
  #   if(is.null(model_switch)) model_switch <- poped.db$downsized.design$model_switch
  #   if(is.null(ni)) ni <- poped.db$downsized.design$ni
  #   if(is.null(xt)) xt <- poped.db$downsized.design$xt
  #   if(is.null(x)) x <- poped.db$downsized.design$x
  #   if(is.null(a)) a <- poped.db$downsized.design$a
  #   if(is.null(groupsize)) groupsize <- poped.db$downsized.design$groupsize
  #   
  if(is.null(model_switch)) model_switch <- poped.db$design$model_switch
  if(is.null(ni)) ni <- poped.db$design$ni
  if(is.null(xt)) xt <- poped.db$design$xt
  if(is.null(x)) x <- poped.db$design$x
  if(is.null(a)) a <- poped.db$design$a
  if(is.null(groupsize)) groupsize <- poped.db$design$groupsize
  
  if(!is.null(fim.calc.type)) poped.db$settings$iFIMCalculationType=fim.calc.type
  if(!is.null(approx.method)) poped.db$settings$iApproximationMethod=approx.method
  if(!is.null(FOCE.num)) poped.db$settings$iFOCENumInd=FOCE.num
  
  if(!is.null(deriv.type)){ 
    poped.db$settings$m1_switch=deriv.type
    poped.db$settings$m2_switch=deriv.type
    poped.db$settings$hle_switch=deriv.type
    poped.db$settings$gradff_switch=deriv.type
    poped.db$settings$gradfg_switch=deriv.type
  }

  output = mftot(model_switch,groupsize,ni,xt,x,a,bpop.val,d_full,sigma_full,docc_full,poped.db)
  FIM <- output$ret
  
  return(FIM)
}
