#'
#'
#' Multiple starts for Regularized Structural Equation Modeling
#'
#'
#' @param model Lavaan output object. This is a model that was previously
#'        run with any of the lavaan main functions: cfa(), lavaan(), sem(),
#'        or growth(). It also can be from the efaUnrotate() function from
#'        the semTools package. Currently, the parts of the model which cannot
#'        be handled in regsem is the use of multiple group models, missing
#'        other than listwise, thresholds from categorical variable models,
#'        the use of additional estimators other than
#'        ML, most notably WLSMV for categorical variables. Note: the model
#'        does not have to actually run (use do.fit=FALSE), converge etc...
#'        regsem() uses the lavaan object as more of a parser and to get
#'        sample covariance matrix.
#' @param max.try number of starts to try before convergence.
#' @param lambda Penalty value. Note: higher values will result in additional
#'        convergence issues.
#' @param LB lower bound vector. Note: This is very important to specify
#'        when using regularization. It greatly increases the chances of
#'        converging.
#' @param UB Upper bound vector
#' @param type Penalty type. Options include "none", "lasso", "ridge",
#'        and "diff_lasso". diff_lasso penalizes the discrepency between
#'        parameter estimates and some pre-specified values. The values
#'        to take the deviation from are specified in diff_par.
#' @param optMethod Solver to use. Recommended options include "nlminb" and
#'        "optimx". Note: for "optimx", the default method is to use nlminb.
#'        This can be changed in subOpt.
#' @param gradFun Gradient function to use. Recommended to use "ram",
#'        which refers to the method specified in von Oertzen & Brick (2014).
#'        The "norm" procedure uses the forward difference method for
#'        calculating the hessian. This is slower and less accurate.
#' @param pars_pen Parameter indicators to penalize. If left NULL, by default,
#'        all parameters in the \emph{A} matrix outside of the intercepts are
#'        penalized when lambda > 0 and type != "none".
#' @param diff_par Parameter values to deviate from. Only used when
#'        type="diff_lasso".
#' @param hessFun Hessian function to use. Recommended to use "ram",
#'        which refers to the method specified in von Oertzen & Brick (2014).
#'        The "norm" procedure uses the forward difference method for
#'        calculating the hessian. This is slower and less accurate.
#' @param verbose Whether to print iteration number.
#' @param warm.start Whether start values are based on previous iteration.
#'        This is not recommended.
#' @param Start2 Provided starting values. Not required
#' @param tol absolute tolerance for convergence.
#' @param max.iter max iterations for optimization.
#' @keywords multiple optim
#' @export
#' @examples
#' \dontrun{
#'
#' library(regsem)
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#' outt = cfa(mod,HS,meanstructure=TRUE)
#'
#' fit1 <- multi_optim(outt,max.try=40,
#'                    lambda=0.1,type="lasso",
#'                    gradFun="ram")
#'
#'
#'# growth model
#'model <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
#'           s =~ 0*t1 + s1*t2 + s2*t3 + 3*t4 '
#'fit <- growth(model, data=Demo.growth)
#'summary(fit)
#'fitmeasures(fit)

#'fit3 <- multi_optim(fit,lambda=0.2,type="ridge",gradFun="none")
#'summary(fit3)
#'}


multi_optim <- function(model,max.try=10,lambda,
                         LB=-Inf,UB=Inf,type,optMethod="nlminb",gradFun="ram",
                         pars_pen=NULL,diff_par=NULL,hessFun="none",
                        verbose=TRUE,warm.start=FALSE,Start2=NULL,
                        tol=1e-6,max.iter=50000){



  if(gradFun=="norm"){
    stop("Only recommended grad function is ram or none at this time")
  }

  if(type=="ridge" & gradFun != "none"){
    warning("At this time, only gradFun=none recommended with ridge penalties")
  }

  if(type=="lasso" & gradFun != "ram"){
    warning("At this time, only gradFun=ram recommended with lasso penalties")
  }


if(warm.start==TRUE){
  stop("warm start not currently functioning")
}
 #   if(warm.start==TRUE){
 #     fit99 <- suppressWarnings(regsem(model,lambda=lambda,type=type,optMethod=optMethod,
  #                                     Start=start.optim,gradFun=gradFun,hessFun=hessFun,max.iter=max.iter,
  #                                     LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par,tol=tol))
 #     start.vals <- as.vector(fit99$coefficients)
#
 #     if(length(start.vals) != length(extractMatrices(model)$parameters)){
 #       start.vals = rep(0.5,length(extractMatrices(model)$parameters)) +
 #         rnorm(length(extractMatrices(model)$parameters),0,0.05)
  #    }

  #  }

  sss = seq(1,max.try)

    mult_run <- function(model,n.try=1,lambda,LB=-Inf,
                         UB=Inf,type,optMethod,
                         gradFun,n.optim,pars_pen,
                         diff_par,hessFun,Start2){
      mtt = matrix(NA,n.try,3)
      mtt[,3] = round(runif(n.try,1,99999))
      start.optim=NULL
      n.optim = 0
      while(n.optim < n.try){
        n.optim=n.optim+1
        set.seed(mtt[n.optim,3])

        if(warm.start==FALSE){
          if(is.null(Start2)){
            start.optim = rep(0.5,length(extractMatrices(model)$parameters)) +
              rnorm(length(extractMatrices(model)$parameters),0,0.05)
          }else{
            start.optim = Start2 +
              rnorm(length(Start2),0,0.05)
          }

        }#else if(warm.start==TRUE){
         # start.optim <- rep(0,length(start.vals))
         # for(i in 1:length(start.vals)){
         #   start.optim[i] <- as.numeric(start.vals[i]) + rnorm(1,0,0.02)
         # }
       # }


        fit1 <- suppressWarnings(try(regsem(model,lambda=lambda,type=type,optMethod=optMethod,
                                            Start=start.optim,gradFun=gradFun,hessFun=hessFun,max.iter=max.iter,
                                            LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par,tol=tol),silent=T))

        if(inherits(fit1, "try-error")) {
          mtt[n.optim,1] = NA
          mtt[n.optim,2] = NA
        }else{
          if(optMethod=="nlminb"){
            mtt[n.optim,1] = fit1$out$objective
            mtt[n.optim,2] = fit1$out$convergence
          }else{
            mtt[n.optim,1] = fit1$out$value
            mtt[n.optim,2] = fit1$out$convcode
          }
        }
      }
      mtt
    }

    iter.optim = 0
    while(iter.optim < max.try){
    iter.optim = iter.optim + 1

    outt = mult_run(model,n.try=1,lambda=lambda,LB,UB,type,
                    optMethod=optMethod,gradFun=gradFun,n.optim=iter.optim,Start2=Start2,
                    pars_pen=pars_pen,diff_par=diff_par,hessFun=hessFun)

    if(verbose==TRUE) print(c(iter.optim,outt[,1],outt[,2]))

      if(all(is.na(outt[,2])==TRUE)){
        return
      }else if(any(na.omit(outt[,2]) == 0)){
        if(any(is.na(outt[,1]) == FALSE & outt[,1] < 9999999999 & outt[,1] > 0)){
        row = which(outt[,1] == min(outt[which(is.na(outt[,1]) == FALSE & outt[,1] > 0 & outt[,2] == 0)]))[1]
        set.seed(outt[row,3])
        if(is.null(Start2)){
          start.optim = rep(0.5,length(extractMatrices(model)$parameters)) +
            rnorm(length(extractMatrices(model)$parameters),0,0.05)
        }else{
          start.optim = Start2 +
            rnorm(length(Start2),0,0.05)
        }

        fit1 <- suppressWarnings(regsem(model,lambda=lambda,type=type,optMethod=optMethod,
                       Start=start.optim,gradFun=gradFun,hessFun=hessFun,max.iter=max.iter,
                       LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par,tol=tol))
        return(fit1)
        break
        }else{
          return
        }
      }else{
        return
      }
    }
   # if(exists("fit1")==FALSE){
      fit1 <- suppressWarnings(regsem(model,lambda=lambda,type=type,optMethod=optMethod,
                     Start="default",gradFun=gradFun,hessFun=hessFun,tol=tol,
                     LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par))

        fit1$convergence <- 99
        return(fit1)
   # }



    #if(fit1$convergence != 0){
      warning("WARNING: Model did not converge! It is recommended to increase max.try")
    #}

}



