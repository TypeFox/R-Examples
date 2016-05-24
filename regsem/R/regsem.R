#'
#'
#' Regularized Structural Equation Modeling
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
#' @param lambda Penalty value. Note: higher values will result in additional
#'        convergence issues. If using values > 0.1, it is recommended to use
#'        mutli_optim() instead. See \code{\link{multi_optim}} for more detail.
#' @param alpha Mixture for elastic net. Not currently working applied.
#' @param type Penalty type. Options include "none", "lasso", "ridge",
#'        and "diff_lasso". diff_lasso penalizes the discrepency between
#'        parameter estimates and some pre-specified values. The values
#'        to take the deviation from are specified in diff_par.
#' @param data Optional dataframe. Only required for missing="fiml" which
#'        is not currently working well.
#' @param optMethod Solver to use. Recommended options include "nlminb" and
#'        "optimx". Note: for "optimx", the default method is to use nlminb.
#'        This can be changed in subOpt.
#' @param gradFun Gradient function to use. Recommended to use "ram",
#'        which refers to the method specified in von Oertzen & Brick (2014).
#'        The "norm" procedure uses the forward difference method for
#'        calculating the hessian. This is slower and less accurate.
#' @param hessFun Hessian function to use. Recommended to use "ram",
#'        which refers to the method specified in von Oertzen & Brick (2014).
#'        The "norm" procedure uses the forward difference method for
#'        calculating the hessian. This is slower and less accurate.
#' @param parallel Logical. Whether to parallelize the processes?
#' @param Start type of starting values to use. Only recommended to use
#'        "default". This sets factor loadings and variances to 0.5.
#'        Start = "lavaan" uses the parameter estimates from the lavaan
#'        model object. This is not recommended as it can increase the
#'        chances in getting stuck at the previous parameter estimates.
#' @param subOpt Type of optimization to use in the optimx package.
#' @param longMod If TRUE, the model is using longitudinal data? This changes
#'        the sample covariance used.
#' @param optNL Type of optimization to use in the NLopt package. Currently
#'        not in use.
#' @param fac.type Using "cfa" or "efa" type of model?
#' @param matrices Function to use for extracting RAM matrices.Only recommended
#'        to use "extractMatrices".
#' @param pars_pen Parameter indicators to penalize. If left NULL, by default,
#'        all parameters in the \emph{A} matrix outside of the intercepts are
#'        penalized when lambda > 0 and type != "none".
#' @param diff_par Parameter values to deviate from. Only used when
#'        type="diff_lasso".
#' @param LB lower bound vector. Note: This is very important to specify
#'        when using regularization. It greatly increases the chances of
#'        converging.
#' @param UB Upper bound vector
#' @param calc Type of calc function to use with means or not. Not recommended
#'        for use.
#' @param tol Absolute tolerance for convergence.
#' @param max.iter Max iterations for optimization.
#' @param missing How to handle missing data. Current options are "listwise"
#'        and "fiml". "fiml" is not currently working well.
#' @return out List of return values from optimization program
#' @return convergence Convergence status.
#' @return par.ret Final parameter estimates
#' @return Imp_Cov Final implied covariance matrix
#' @return grad Final gradient.
#' @return KKT1 Were final gradient values close enough to 0.
#' @return KKT2 Was the final Hessian positive definite.
#' @return df Final degrees of freedom. Note that df changes with lasso
#'         penalties.
#' @return npar Final number of free parameters. Note that this can change
#'         with lasso penalties.
#' @return SampCov Sample covariance matrix.
#' @return fit Final F_ml fit. Note this is the final parameter estimates
#'         evaluated with the F_ml fit function.
#' @return coefficients Final parameter estimates
#' @return nvar Number of variables.
#' @return N sample size.
#' @return nfac Number of factors
#' @return baseline.chisq Baseline chi-square.
#' @return baseline.df Baseline degrees of freedom.
#' @keywords optim calc
#' @useDynLib regsem
#' @import Rcpp
#' @import lavaan
#' @importFrom stats cov na.omit nlminb pchisq rnorm runif sd uniroot var
#' @export
#' @examples
#' library(lavaan)
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#' # Recommended to specify meanstructure in lavaan
#' outt = cfa(mod,HS,meanstructure=TRUE)
#'
#' fit1 <- regsem(outt,lambda=0.1,type="lasso",gradFun="ram")





regsem = function(model,lambda=0,alpha=0,type="none",data=NULL,optMethod="nlminb",
                 gradFun="ram",hessFun="none",parallel="no",Start="default",
                 subOpt="nlminb",longMod=F,
                 optNL="NLOPT_LN_NEWUOA_BOUND",fac.type="cfa",
                 matrices="extractMatrices",
                 pars_pen=NULL,
                 diff_par=NULL,
                 LB=-Inf,
                 UB=Inf,
                 calc="normal",
                 tol=1e-6,
                 max.iter=150,
                 missing="listwise"){


    if(gradFun=="norm"){
      stop("Only recommended grad function is ram or none at this time")
    }

  if(type=="ridge" & gradFun != "none"){
    warning("At this time, only gradFun=none recommended with ridge penalties")
  }

  if(type=="lasso" & gradFun != "ram"){
    warning("At this time, only gradFun=ram recommended with lasso penalties")
  }

    parL = parTable(model)[,"label"]
    if(sum(duplicated(parL[parL != ""])) > 0){
      stop("regsem currently does not allow equality constraints")
    }

    nvar = model@pta$nvar[[1]][1]
    nfac = model@pta$nfac[[1]][1]

    if(missing=="listwise"){
      calc_fit = "cov"
      nobs = model@SampleStats@nobs[[1]][1]
      SampMean <- model@SampleStats@mean[][[1]]

      if(length(model@ParTable$op[model@ParTable$op == "~1"]) > 0){
        SampCov1 <- model@SampleStats@cov[][[1]]
        SampCov2 <- SampCov1 + SampMean%*%t(SampMean)
        # try changing size of SampCov
        SampCov3 = cbind(SampCov2,SampMean)
        SampCov = rbind(SampCov3,append(SampMean,1))
      }else if(length(model@ParTable$op[model@ParTable$op == "~1"]) == 0){
        SampCov <- model@SampleStats@cov[][[1]]
      }

      #for grad ram with mean
      SampCov22 <- model@SampleStats@cov[][[1]]# + SampMean %*% t(SampMean)

    }else if(missing=="fiml"){
      #stop("FIML is currently not supported at this time")
      calc_fit = "ind"
      if(is.null(data)==TRUE){
        stop("Dataset needs to be provided for missing==fiml")
      }


      if(length(model@ParTable$op[model@ParTable$op == "~1"]) == 0){
        stop("meanstructure needs to be equal to TRUE for FIML")
      }

    }
    #SampCov <- fitted(model)$cov

    #SampMean <- rep(0,nvar)

    type2 = 0
    if(type=="lasso"){
      type2 = 1
    }else if(type=="ridge"){
      type2=2
    }else if(type=="diff_lasso"){
      type2=3
    }




    #nUniq = nvar
    #nFacCov
    df = model@Fit@test[[1]]$df
    npar = model@Fit@npar
    nload = length(model@ParTable$op[model@ParTable$op == "=~"])


if(matrices == "semPlot"){
if(fac.type=="cfa"){
    # get matrices using semPlot package
    out = semPlot::modelMatrices(model)

    A = out$A[[1]]$par
    S = out$S[[1]]$par
    F = out$F[[1]]$est
    A_fixed = out$A[[1]]$fixed
    #A.std = out$A[[1]]$std # try not to use
    A_est = out$A[[1]]$est
    S_fixed = out$S[[1]]$fixed
    #S.std = out$S[[1]]$std
    S_est = out$S[[1]]$est
    I <- diag(nrow(A))

} else if(fac.type=="efa"){
      out = semPlot::modelMatrices(model)
      # if efaUnrotate minus S number by number of variables for intercepts
      A = out$A[[1]]$par
      S = out$S[[1]]$par
      S[S != 0] <- S[S != 0] - nvar

      F = out$F[[1]]$est
      A.fixed = out$A[[1]]$fixed
      A.std = out$A[[1]]$std # try not to use
      A.est = out$A[[1]]$est
      S.fixed = out$S[[1]]$fixed
      S.std = out$S[[1]]$std
      S.est = out$S[[1]]$est
      I <- diag(nrow(A))
}
}else if(matrices == "extractMatrices"){
  list <- extractMatrices(model)
  A <- list$A
  A_est <- list$A_est
  A_fixed <- list$A_fixed
  S <- list$S
  S_est <- list$S_est
  S_fixed <- list$S_fixed
  F <- list$F
  I <- diag(nrow(A))
}


    if(is.null(pars_pen) == TRUE){
      if(any(colnames(A) == "1")){
        IntCol = which(colnames(A) == "1")
        A_minusInt = A[,-IntCol]
        A_pen = A_minusInt != 0
        pars_pen = A_minusInt[A_pen]
      }else{
        A_pen = A != 0
        pars_pen = A[A_pen]
      }
    }

    # set fixed parameters
    ############## probably should be set after estimates are given as parameters fixed to
    ############## 1 will be thought to be start parameter 1
    #A[A.fixed ==T] <- A.est[A.fixed==T]
    #S[S.fixed ==T] <- S.est[S.fixed==T]

   if(class(Start)=="numeric"){
      start=Start
   }else if(class(Start) != "numeric"){
     if(Start=="lavaan"){
       # get starting values
       lambda.start = lavaan::inspect(model,"start")$lambda
       psi.start = lavaan::inspect(model,"start")$psi  # tricky to get ordering
       theta.start = diag(lavaan::inspect(model,"start")$theta) # only diagonal elements of theta
       nu.start = lavaan::inspect(model,"start")$nu
       alpha.start = lavaan::inspect(model,"start")$alpha
       # put into vector
       #par.start = as.vector(c(lambda.start[lambda.start != 0],psi.start[!upper.tri(psi.start)],theta.start),mode="numeric")
       par.start = as.vector(c(lambda.start,theta.start,psi.start[!upper.tri(psi.start)]),mode="numeric")
       # assign which values in parTable have a number assigned !=0 in free column
       free <- lavaan::parTable(model)$free
       # pull only free parameters
       start = par.start[free > 0]
       st <- is.na(start)
       start <- start[!st]

     } else if(Start == "default"){
       nstart <- max(max(A),max(S))
       start <- rep(0.5,nstart)

     }
   }else if(Start=="prev"){
     nstart <- max(max(A),max(S))
     start.seq <- seq(1:nstart)
     start <- rep(0.5,nstart)

     for(i in 1:nstart){
       if(sum(A == i) > 0) {
         start[i] = A.est[A == i]
       }
       else if(sum(S== i)>0){
         start[i] = S.est[S==i][1]
       }
     }
   }


  if(calc == "normal"){
    calc = function(start){
         mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
         #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
         pen_vec = c(mult$A_est22[A %in% pars_pen],mult$S_est22[S %in% pars_pen])
         if(type=="diff_lasso"){
           pen_diff = pen_vec - diff_par
         }else{
           pen_diff=0
         }
         if(calc_fit=="cov"){
           #fit = fit_fun(ImpCov=mult$ImpCov,SampCov,Areg=mult$A_est22,lambda,alpha,type,pen_vec)
           fit = rcpp_fit_fun(ImpCov=mult$ImpCov,SampCov,type2,lambda,pen_vec,pen_diff)
           fit
         }else if(calc_fit=="ind"){
           #stop("Not currently supported")
           fit = fiml_calc(ImpCov=mult$ImpCov,data=data,
                           Areg=mult$A_est22,lambda,alpha,type,pen_vec,nvar)
         }

    }
  }else if(calc == "calc2"){
    calc = function(start){
      #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
      fit = ram_calc(par=start,SampCov22,A,S,F,SampMean)
      like = fit$lik
      like
    }
  }else if(calc == "boot"){
    calc = function(start){
      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
      pen_vec = c(mult$A_est22[A %in% pars_pen],mult$S_est22[S %in% pars_pen])
      if(type=="diff_lasso"){
        pen_diff = pen_vec - diff_par
      }else{
        pen_diff=0
      }
      n.boot=100
      fits = rep(NA,n.boot)
      # boot part
      for(i in 1:n.boot){
        dat1 = model@Data@X[[1]]
        ids <- sample(nrow(dat1),nrow(dat1),replace=TRUE)
        new.dat <- dat1[ids,]
        SampCov.boot <- cov(new.dat)
        #fits[i] = rcpp_fit_fun(ImpCov=mult$ImpCov,SampCov.boot,type2,lambda,pen_vec,pen_diff)
        fits[i] = fit_fun(ImpCov=mult$ImpCov,SampCov.boot,lambda,alpha,type,pen_vec)
      }
      mean(fits)
    }
  }
#------------------------ only works for CFA models --------------------------------
#  grad = function(start){
#    mult = RAM_multSimp(A,A.fixed,S,S.fixed,F,start,nfac,nvar)
#    ret = gradient(ExpCov=mult$ExpCov,cov,A=mult$A,S=mult$S,
#                   start,lambda,alpha,type,nvar,nload,nfac,nUniq,nFacCov)
#    ret
# }

  if(gradFun=="norm"){
    grad = function(start){

      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
                #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
                ret = grad_fun(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                               Sreg=mult$S_est22,A,A_fixed,A_est,S,S_fixed,S_est,
                               F,lambda,alpha,type,pars_pen)
               ret
  }
} else if(gradFun=="ram"){
    grad = function(start){

      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
     #   ret = grad_ram(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
     #                  Sreg=mult$S_est22,A,S,
      #                  F,lambda,type,pars_pen,diff_par)
      #pen_vec = c(mult$A_est22[A %in% pars_pen],mult$S_est22[S %in% pars_pen])
       ret = rcpp_grad_ram(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                       Sreg=mult$S_est22,A,S,
                         F,lambda,type2=type2,pars_pen,diff_par=0)
      ret
    }

}else if(gradFun=="ram_mean"){
  grad = function(start){

    mult = ram_calc(par=start,SampCov22,A,S,F,SampMean)
    ret = grad_ram_wMean(par=start,ImpCov=mult$ImpCov,SampCov22,Areg = mult$A2,
                    Sreg=mult$S2,A=mult$A.pars,S=mult$S.pars,
                    F=mult$F,SampMean,lambda,type,m=mult$m,mu=mult$mu,m.pars=mult$m.pars)
    ret
  }
} else if(gradFun=="numDeriv"){
    grad = function(start){

    #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
    ret = numDeriv::grad(calc,x=start)
    ret
  }
}


  if(hessFun=="norm"){

    if(parallel=="no"){

    hess = function(start){

        mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
        #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
        retH = hessian(par=start,ImpCov=mult$ImpCov,SampCov,A,A_fixed,A_est,
                 S,S_fixed,S_est,F)
        retH
    }
  } else if(parallel=="yes"){

    hess = function(start){
    mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
    #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
    retH = hessian_parallel(par=start,ImpCov=mult$ImpCov,A,A.fixed,A.est,
                   S,S.fixed,S.est,F)
    retH
  }
}
}  else if(hessFun=="ram"){

    if(parallel=="no"){
    hess = function(start){

      mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      #mult = RAMmult(par=start,A,S,F,A_fixed,A_est,S_fixed,S_est)
      ret = hess_ram(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A_est22,
                      Sreg=mult$S_est22,A,S,F)
      ret
    }
  } else if(parallel=="yes"){
      hess = function(start){
        mult = rcpp_RAMmult(par=start,A,S,S_fixed,A_fixed,A_est,S_est,F,I)
        #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
        ret = hess_ramParallel(par=start,ImpCov=mult$ImpCov,SampCov,Areg = mult$A.est22,
                        Sreg=mult$S.est22,A,S,F)
        ret
      }

    }
} else if(hessFun=="numDeriv"){
  hess = function(start){

    #mult = RAMmult(par=start,A,S,F,A.fixed,A.est,S.fixed,S.est)
    ret = numDeriv::hessian(calc,x=start)
    ret
  }
}else{
  hess = NULL
}

    res <- list()

if(optMethod=="nlminb"){
    if(gradFun=="norm"){
      if(hessFun=="norm"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(abs(diag(S)-max(A)),rep(-10,max(S)-max(diag(S))))
        #UB = c(100,100,100,1)
        out <- nlminb(start,calc,grad,hess,lower=LB,upper=UB,control=list(eval.max=max.iter,
                                                                 iter.max=max.iter))
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
       #LB = c(rep(-6,max(A)),rep(1e-6,rep(-10,max(S)-max(diag(S))),max(diag(S))-max(A)))
        out <- nlminb(start,calc,grad,lower=LB,upper=UB,control=list(eval.max=max.iter,
                                                                     iter.max=max.iter))
        res$out <- out
        res$convergence = out$convergence
        #res$optim_fit <- out$objective
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }else if(gradFun=="ram"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
       out <- nlminb(start,calc,grad,hess,lower=LB,upper=UB,control=list(eval.max=max.iter,
                                                                         iter.max=max.iter))
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
       out <- nlminb(start,calc,grad,lower=LB,upper=UB,
                     control=list(eval.max=max.iter,
                     iter.max=max.iter,step.min=0.0000001)) #,x.tol=1.5e-6
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }else if(gradFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,lower=LB,upper=UB,control=list(eval.max=max.iter,
                                                                 iter.max=max.iter,
                                                                step.min=0.0000001))
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
    }else if(gradFun=="ram_mean"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,grad,hess,lower=LB,control=list(eval.max=max.iter,
                                                                 iter.max=max.iter))
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc2,grad,lower=LB,upper=UB,eval.max=max.iter,
                      iter.max=max.iter)
        res$out <- out
        res$convergence = out$convergence
        #res$optim_fit <- out$objective
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }else if(gradFun=="none"){
      #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
      out <- nlminb(start,calc,lower=LB,upper=UB,eval.max=max.iter,
                    iter.max=max.iter)
      res$out <- out
      #res$optim_fit <- out$objective
      res$convergence = out$convergence
      par.ret <- out$par
      #res$iterations <- out$iterations
    }else if(gradFun=="numDeriv"){
      if(hessFun=="numDeriv"){
        warning("numDeriv does not seem to be accurate at this time")
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,grad,hess,lower=LB,upper=UB,eval.max=max.iter,
                      iter.max=max.iter)
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }else if(hessFun=="none"){
        warning("numDeriv does not seem to be accurate at this time")
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- nlminb(start,calc,lower=LB,upper=UB,gradient=grad,eval.max=max.iter,
                      iter.max=max.iter)
        res$out <- out
        #res$optim_fit <- out$objective
        res$convergence = out$convergence
        par.ret <- out$par
        #res$iterations <- out$iterations
      }
    }
}else if(optMethod=="optimx"){
    if(gradFun=="norm"){
      if(hessFun=="norm"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,hess,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals

        res$convergence = out$convcode
        par.ret <- coef(out)
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,lower=LB,upper=UB,method=subOpt,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals
        res$convergence = out$convcode
        par.ret <- coef(out)
      }
    }else if(gradFun=="ram_mean"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,hess,lower=LB,upper=UB,method=subOpt,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals
        res$convergence = out$convcode
        par.ret <- coef(out)
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE))
        res$out <- out
        res$convergence = out$convcode
        par.ret <- coef(out)
       }
    }else if(gradFun=="ram"){
      if(hessFun=="ram"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,hess,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE,itnmax = max.iter))
        res$out <- out
        #res$optim_fit <- out$value
        res$convergence = out$convcode
        par.ret <- coef(out)
      }else if(hessFun=="none"){
        #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-max(A)),rep(-10,max(S)-max(diag(S))))
        out <- optimx::optimx(start,calc,grad,method=subOpt,lower=LB,upper=UB,control=list(starttests=FALSE))
        res$out <- out
        #res$iterations <- out$fevals
        #res$optim_fit <- out$value
        res$convergence = out$convcode
        par.ret <- coef(out)
      }
    }else if(gradFun=="none"){
      #LB = c(rep(-6,max(A)),rep(1e-6,max(diag(S))-S[1,1]+1),rep(-10,max(S)-max(diag(S))))
      out <- optimx::optimx(start,calc,method=subOpt,itnmax=1500,lower=LB,upper=UB,
                    hessian=T,control=list(maxit=1500,starttests=FALSE,all.methods=TRUE,abstol=1e-12))
      res$out <- out
      #res$iterations <- out$fevals
      res$convergence = out$convcode
      #pars <- coef(out)
      #res$pars <- pars
      par.ret <- coef(out)
    }
}else if(optMethod=="rsolnp"){
       # if(UB == Inf) UB=NULL
       # if(LB == -Inf) LB=NULL
        suppressWarnings(out <- Rsolnp::solnp(start,calc,LB=LB,UB=UB,
                                              control=list(trace=0,tol=1e-16)))#tol=1e-16
        #out <- optim(par=start,fn=calc,gr=grad)
        res$out <- out
        #res$iterations <- out$nfunevals
        res$convergence = out$convergence
        par.ret <- out$pars

}else if(optMethod=="rgenoud"){
  dom = matrix(c(LB,UB),nrow=length(start),2)
  suppressWarnings(out <- rgenoud::genoud(calc,starting.values=start,Domains=dom,
                                      nvars=length(start),print.level=0,gr=grad,
                                      boundary.enforcement=2,
                                      wait.generations=3))
  res$out <- out
  res$optim_fit <- out$value
  res$convergence = 0
  res$par.ret <- out$par
}else if(optMethod=="GA"){
  calc2 = function(start){
    10-calc(start)
  }
  out = GA::ga("real-valued", fitness = calc2, nBits = length(start),
               min=LB,max=UB,monitor=FALSE,pcrossover=0.9,pmutation=0.3,
               maxiter=10000)
  res$out <- summary(out)
  res$optim_fit <- 10 - summary(out)$fitness
  res$convergence = 0
  res$par.ret <- summary(out)$solution
}else if(optMethod=="optim_rj"){
  out = optim_rj(start=start,func=calc,grad=grad,hess=hess,tol=tol,max.iter=max.iter)
  res$out <- out
  res$optim_fit <- out$value
  res$convergence = out$convergence
  res$par.ret <- out$pars
  res$iterations <- out$iterations
}





   # imp_cov = RAMmult(res$out$par,A,S,F,A.fixed,A.est,S.fixed,S.est) [[1]]
   # res$imp_cov <- imp_cov
    #res$df <- df
    #res$nobs <- nobs
    #res$nload <- nload

    #hess <- ret_hess(res$out$par,A,S,F,A_fixed,A_est,S_fixed,S_est)
    #res$hess <- hess
    #res$info <- ginv((nobs/2)*hess)
    #res$A <- A
    #res$S <- S
    #res$A.est <- A_est
    #res$A.fixed <- A_fixed
    #res$S_est <- S_est
    #res$S.fixed <- S_fixed
   # res$F <- F


    pars.df <- data.frame(matrix(NA,1,max(max(A),max(S))))
    pars.df[1,] <- par.ret

    if(any(pars.df[diag(S[diag(S) != 0])] < 0)){
      warning("Some Variances are Negative!")
    }


    #res$ftt = rcpp_RAMmult(par=as.numeric(pars.df),A,S,S_fixed,A_fixed,A_est,S_est,F,I)
      # get Implied Covariance matrix

    Imp_Cov1 <- rcpp_RAMmult(par=as.numeric(pars.df),A,S,S_fixed,A_fixed,A_est,S_est,F,I)$ImpCov
    #Imp_Cov <- RAMmult(par=as.numeric(pars.df),A,S,F,A_fixed,A_est,S_fixed,S_est)$ImpCov



    if(length(model@ParTable$op[model@ParTable$op == "~1"]) > 0 & missing=="listwise"){
      Imp_Cov = Imp_Cov1[1:(nrow(Imp_Cov1)-1),1:(ncol(Imp_Cov1)-1)] - SampMean %*% t(SampMean)
    }else{
      Imp_Cov = Imp_Cov1
    }

    res$Imp_Cov <- Imp_Cov

    #res$grad <- grad(as.numeric(pars.df))
    #### KKT conditions #####
    if(gradFun=="none"){
      res$KKT1 = "grad not specified"
    }else{
      res$grad <- grad(as.numeric(pars.df))
      kk = try(all(grad(as.numeric(pars.df)) < 0.001))
      if(inherits(kk, "try-error")){
        res$KKT1 = "error"
      }else{
        if(kk == TRUE){
          res$KKT1 = TRUE
        }else if(kk < 0.001){
          res$KKT1 = FALSE
        }else{
          res$KKT1 = NA
        }
      }
    }


    if(hessFun=="none"){
      res$KKT2 = "hess not specified"
    }else{
      hess.mat = hess(as.numeric(pars.df))
      eig = eigen(hess.mat)$values
      hess.eigs = try(all(eig) > 1e-6)
      if(inherits(hess.eigs, "try-error")){
        res$KKT2 = "error"
      }else{
        if(hess.eigs == TRUE){
          res$KKT2 = TRUE
        }else if(hess.eigs == FALSE){
          res$KKT2 = FALSE
        }else{
          res$KKT2 = NA
        }
      }
    }




   # rettt = rcpp_RAMmult(par=as.numeric(pars.df),A,S,S_fixed,A_fixed,A_est,S_est,F,I)
   # A_new2 <- rettt$A_est22
   # S_new2 <- rettt$S_est22
    #A_new <- A
    #S_new <- S
    #I = diag(ncol(A))

    #for(i in 1:max(A)) A_new[A_new==i] = pars.df[i]
    #A_new2 = matrix(unlist(A_new),nrow(A),ncol(A))

    #for(i in min(S[S>0]):max(S)) S_new[S_new==i] = pars.df[i]
    #S_new2 = matrix(unlist(S_new),nrow(S),ncol(S))
    #res$retttt = calc(as.numeric(pars.df))
    #res$Imp_Cov = F %*% solve(I-A_new2) %*% S_new2 %*% t(solve(I-A_new2)) %*% t(F)
    #res$fitt =  calc(as.numeric(pars.df))

    for(i in 1:ncol(pars.df)){
      if(any(A == i)){
        pos = which(A == i,arr.ind=T)
        one = colnames(A)[pos[2]]
        two = rownames(A)[pos[1]]
        names(pars.df)[i] = paste(one,"->",two)
      }else if(any(S==i)){
        pos = which(S == i,arr.ind=T)
        one = colnames(S)[pos[2]]
        two = rownames(S)[pos[1]]
        names(pars.df)[i] = paste(one,"~~",two)
      }
    }


    if(type=="none" | lambda==0){
      res$df = df
      res$npar = npar
    }else if(type=="lasso" | type=="ridge"){
      #A_estim = A != 0
      #pars = A_est[A_estim]
      pars_sum = pars.df[pars_pen]
      pars_l2 = sqrt(pars_sum**2)
      res$df = df + sum(pars_l2 < 0.001)
      res$npar = npar - sum(pars_l2 < 0.001)

    }
    if(missing == "listwise"){
     # SampCov <- model@SampleStats@cov[][[1]]
    #  res$SampCov = SampCov
      #res$fit = 0.5*(log(det(Imp_Cov1)) + trace(SampCov %*% solve(Imp_Cov1)) -
       #              log(det(SampCov))  - nvar)
      res$fit = rcpp_fit_fun(Imp_Cov1, SampCov,type2=0,lambda=0,pen_vec=0,pen_diff=0)
    }

    SampCov <- model@SampleStats@cov[][[1]]
    res$SampCov = SampCov

    res$coefficients <- round(pars.df,3)
    res$nvar = nvar
    res$N = nobs
    res$nfac = nfac

    if(model@Fit@converged == FALSE){
      res$baseline.chisq = NA
      res$baseline.df = NA
    }else{
      res$baseline.chisq = fitMeasures(model)["baseline.chisq"]
      res$baseline.df = fitMeasures(model)["baseline.df"]
    }

    #res$grad <- grad(res$par.ret)

    #res$hess <- hess(res$par.ret)


    if(res$convergence != 0){
      warning("WARNING: Model did not converge! It is recommended to try multi_optim()")
    }

    res$call <- match.call()
    class(res) <- "regsem"
    return(res)
}
