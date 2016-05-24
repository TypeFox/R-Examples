#'
#'
#' The main function that ties together and runs the models.
#' @param model lavaan output object.
#' @param n.lambda number of penalization values to test.
#' @param mult.start Logical. Whether to use multi_optim() (TRUE) or
#'         regsem() (FALSE).
#' @param multi.iter number of random starts for multi_optim
#' @param jump Amount to increase penalization each iteration.
#' @param type penalty type.
#' @param fit.ret Fit indices to return.
#' @param fit.ret2 Return fits from just train sample?
#' @param data Optional dataframe. Only required for missing="fiml".
#' @param optMethod solver to use.
#' @param gradFun gradient function to use.
#' @param hessFun hessian function to use.
#' @param parallel whether to parallelize the processes?
#' @param Start type of starting values to use.
#' @param subOpt type of optimization to use in the optimx package.
#' @param longMod longitudinal model?
#' @param optNL type of optimization to use in the NLopt package.
#' @param fac.type using cfa or efa type of model.
#' @param matrices function to use for extracting RAM matrices.
#' @param pars_pen parameter indicators to penalize.
#' @param diff_par parameter values to deviate from.
#' @param LB lower bound vector.
#' @param UB upper bound vector
#' @param calc type of calc function to use with means or not.
#' @param tol absolute tolerance for convergence.
#' @param max.iter max iterations for optimization.
#' @param missing How to handle missing data. Current options are "listwise"
#'        and "fiml".
#' @param ... Any additional arguments to pass to regsem() or multi_optim().
#' @keywords optim calc
#' @export
#' @examples
#' \dontrun{
#' library(lavaan)
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#' outt = cfa(mod,HS)
#'
#' cv.out = cv_regsem(outt,type="ridge",gradFun="none",n.lambda=100)
#'}



cv_regsem = function(model,
                     n.lambda=100,
                     mult.start=TRUE,
                     multi.iter=100,
                     jump=0.002,
                     type="none",
                     fit.ret=c("rmsea","BIC"),
                     fit.ret2 = c("train","test"),
                     data=NULL,
                     optMethod="nlminb",
                    gradFun="ram",
                    hessFun="none",
                    parallel="no",
                    Start="default",
                    subOpt="nlminb",
                    longMod=F,
                    optNL="NLOPT_LN_NEWUOA_BOUND",
                    fac.type="cfa",
                    matrices="extractMatrices",
                    pars_pen=NULL,
                    diff_par=NULL,
                    LB=-Inf,
                    UB=Inf,
                    calc="normal",
                    tol=1e-10,
                    max.iter=50000,
                    missing="listwise",
                    ...){





par.matrix <- matrix(0,n.lambda,model@Fit@npar)
fits <- matrix(NA,n.lambda,length(fit.ret)+2)
SHRINK = 0
count = 0
counts=n.lambda
#res2 <- data.frame(matrix(NA,counts,3))
#coefs = rep(1,14)

while(count < counts){

  count = count + 1
  SHRINK <- jump*(count-1) # 0.01 works well & 0.007 as well with 150 iterations

if(mult.start==FALSE){
  out <- regsem(model=model,lambda=SHRINK,type=type,data=data,
                   optMethod=optMethod,
                   gradFun=gradFun,hessFun=hessFun,
                   parallel=parallel,Start=Start,
                   subOpt=subOpt,
                   longMod=longMod,
                   optNL=optNL,
                   fac.type=fac.type,
                   matrices=matrices,
                   pars_pen=pars_pen,
                   diff_par=diff_par,
                   LB=LB,
                   UB=UB,
                   calc=calc,
                   tol=tol,
                   max.iter=max.iter,
                   missing=missing)


  }else if(mult.start==TRUE){
   out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                      LB=LB,UB=UB,type=type,optMethod=optMethod,
                      gradFun=gradFun,hessFun=hessFun,
                      pars_pen=pars_pen,diff_par=NULL)
  }


  #if(any(fit.ret2 == "test")==TRUE){
  #  fits[[count]]$test = NA #fit_indices(out,CV=TRUE)[fit.ret]
  #}else
  if(any(fit.ret2 == "train")==TRUE){
    fitt = try(fit_indices(out,CV=FALSE)$fits[fit.ret],silent=T)
    if(inherits(fitt, "try-error")) {
      fits[count,3:ncol(fits)] = rep(NA,ncol(fits)-2)
    }else{
      fits[count,3:ncol(fits)] = fitt
    }

  }
  fits[count,1] <- SHRINK
  fits[count,2] <- out$out$convergence

  if(is.null(out$coefficients)==TRUE){
    break
  }
  par.matrix[count,] = as.matrix(out$coefficients)

}
#fits = fit_indices(out,CV=FALSE)
colnames(par.matrix) = names(out$coefficients)
colnames(fits) <- c("lambda","conv",fit.ret)
ret <- list(par.matrix,fits)

}
