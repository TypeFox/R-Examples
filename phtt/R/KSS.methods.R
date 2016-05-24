## Methods ========================================================================================

KSS <- function(formula,
                additive.effects = c("none", "individual", "time", "twoways"),
                consult.dim.crit = FALSE,
                d.max            = NULL,
                sig2.hat         = NULL,
                factor.dim       = NULL,
                level            = 0.01,
                spar             = NULL,
                CV               = FALSE,
                convergence      = 1e-6,#.Machine$double.eps^0.25,
                restrict.mode    = c("restrict.factors","restrict.loadings"),
                ...){
  UseMethod("KSS")
}

print.KSS <- function(x,...){
   cat("Call:\n")
   print(x$call)

   cat("\nCoeff(s) of the Observed Regressor(s) :\n\n")
   slope.para <- x$slope.para
    if(x$is.intercept){
      inter <- matrix(x$Intercept, 1, 1)
      colnames(inter) <- ""
      rownames(inter) <- "(Intercept)"
      slope.para <- rbind(signif(inter,digits=3), signif(slope.para,digits=3))
    }
   cat(signif(t(slope.para),5))
   cat("\n\nAdditive Effects Type: ", as.name(x$additive.effects)," \n")
   cat("\nDimension of the Unobserved Factors:", x$used.dim," \n")
 }
  
coef.KSS <- function(object,...){
    if(object$is.intercept)
    Intercept <- object$Intercept
    else Intercept <- NULL
    
    Slope.Coef <- object$slope.para
    
    if(object$additive.effects== "individual"| object$additive.effects== "twoways")
    Add.Ind.Eff <- object$Add.Ind.Eff
    else Add.Ind.Eff <- NULL

    if(object$additive.effects== "time"| object$additive.effects== "twoways")
    Add.Tim.Eff <- object$Add.Tim.Eff
    else Add.Tim.Eff <- NULL
    
    Common.factors <- object$unob.factors
    
    Ind.loadings.param <- object$ind.loadings
    
    Time.varying.ind.eff <- object$unob.fact.stru
    
    Factor.Dim <- object$used.dim
    if(Factor.Dim>0){
      Var.shares.of.loadings.param      <- numeric(Factor.Dim)
      Total.var.loadings.param          <- sum(apply(Ind.loadings.param,2,var))
      for(i in 1:Factor.Dim){
        Var.shares.of.loadings.param[i] <- round(var(c(Ind.loadings.param[,i]))/Total.var.loadings.param,
                                                 digits=4)*100
      }
    }else{
      Var.shares.of.loadings.param      <- NULL
    }
    coef.list <- list(
        Intercept                    = Intercept,
        Slope.Coef                   = Slope.Coef,
        Add.Ind.Eff                  = Add.Ind.Eff, 
        Add.Tim.Eff                  = Add.Tim.Eff, 
        Common.factors               = Common.factors, 
        Ind.loadings.param           = Ind.loadings.param,
        Var.shares.of.loadings.param = Var.shares.of.loadings.param,
        Time.varying.ind.eff         = Time.varying.ind.eff,
        Factor.Dim                   = Factor.Dim)  
        
    return(coef.list)
}

residuals.KSS <- resid.KSS <- function(object,...){
  Residual.mat <- object$residuals
  return(Residual.mat)
}

summary.KSS <- function(object,...){
  ## Residuals:
  Res.outpt <- round((summary(as.vector(object$residuals))), digits=2)[-4]
  names(Res.outpt) <- c("Min", "1Q", "Median", "3Q", "Max")
  yy <- sum(diag(crossprod(object$orig.Y - mean(object$orig.Y))))
  ee <- sum(diag(crossprod(object$residuals)))
  R2  <- 1 - ee/yy
  ## R2a <- 1 - (ee/object$degrees.of.freedom)/(yy - 1)
  
  ## Add-Effect-Type:
  eff              <- matrix(object$additive.effects)
  colnames(eff)    <- ""
  rownames(eff)    <- ""
  
  ## Inference for Coefficients: ============================================================
  if(object$is.intercept){
    Intercept.se   <- sqrt(object$Intercept.V)
    beta.se        <- sqrt(diag(object$beta.V))
    Intercept.zval <- coef(object)$Intercept  / Intercept.se
    beta.zval      <- coef(object)$Slope.Coef / beta.se
    TAB  <- cbind("Estimate" = c(signif(as.numeric(coef(object)$Intercept),             digits=3),
                                 signif(as.numeric(coef(object)$Slope.Coef),            digits=3)),
                  "StdErr"   = c(signif(Intercept.se, digits=3),
                                 signif(beta.se,      digits=3)),
                  "z.value"  = c(signif(as.numeric(Intercept.zval),                digits=3),
                                 signif(as.numeric(beta.zval),                     digits=3)),
                  "Pr(>z)"   = c(signif(as.numeric(2*pnorm(-abs(Intercept.zval))), digits=3),
                                 signif(as.numeric(2*pnorm(-abs(beta.zval))),      digits=3))
                  )  
    rownames(TAB)  <- c("(Intercept)", object$names[2:length(object$names)])
  }else{
    beta.se        <- sqrt(diag(object$beta.V))
    beta.zval      <- coef(object)$Slope.Coef / beta.se
    TAB  <- cbind("Estimate" = signif(as.numeric(coef(object)$Slope.Coef),              digits=3),
                  "StdErr"   = signif(beta.se,                                     digits=3),
                  "z.value"  = signif(as.numeric(beta.zval),                       digits=3),
                  "Pr(>z)"   = signif(as.numeric(2*pnorm(-abs(beta.zval))),        digits=3)
                  )  
    rownames(TAB)  <- object$names[2:length(object$names)]
  }
  
  ## Result: ==============================================================================
  result        <- list(Res.outpt    = Res.outpt,
                        coefficients = TAB,
                        R2           = R2,
                        KSS.obj      = object)                
  class(result) <- "summary.KSS"
  result
}


print.summary.KSS <- function(x, ...){
  ## Call
  cat("Call:\n")
  print(x$KSS.obj$call)
  ## Residuals:
  cat("\nResiduals:\n")
  print(x$Res.outpt)
  cat("\n")
  ## Beta-Coeffs
  cat("\n Slope-Coefficients:\n")
  printCoefmat(x$coefficients)
  
  cat("\nAdditive Effects Type: ",                   as.name(x$KSS.obj$additive.effects)," \n")
  cat("\nUsed Dimension of the Unobserved Factors:", x$KSS.obj$used.dim, " \n")
#  cat("\nOptimized Factor Dimension:              ", x$KSS.obj$optimal.dim," \n") 
  cat("\nResidual standard error:",             signif(x$KSS.obj$sig2.hat, digits=3), "on", 
                                                x$KSS.obj$degrees.of.freedom, "degrees of freedom \n")
  cat("R-squared:",                    signif(x$R2,digits=3),"\n")
#  cat("Adjusted R-squared:",                    signif(x$R2,digits=3),"\n")
}


plot.summary.KSS <- function(x,...){
  if(is.null(x$KSS.obj$unob.factors) & x$KSS.obj$additive.effects=="none"){
    stop("Neither an estimated factor structure nor additive effects to plot.")
  }
  if(!is.null(x$KSS.obj$unob.factors)){
    if(x$KSS.obj$additive.effects=="none"){
      par(mfrow=c(1,2))
      matplot(x$KSS.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$KSS.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(x$KSS.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$KSS.obj$additive.effects=="time"){
      par(mfrow=c(1,3))
      plot.ts(x$KSS.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(x$KSS.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$KSS.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(x$KSS.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$KSS.obj$additive.effects=="twoways"){
      par(mfrow=c(1,4))
      plot.ts(x$KSS.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(x$KSS.obj$Add.Ind.Eff,each=x$KSS.obj$dat.dim[1]),
                     nrow=x$KSS.obj$dat.dim[1],ncol=x$KSS.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(x$KSS.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$KSS.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(x$KSS.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$KSS.obj$additive.effects=="individual"){
      par(mfrow=c(1,3))
      matplot(matrix(rep(x$KSS.obj$Add.Ind.Eff,each=x$KSS.obj$dat.dim[1]),
                     nrow=x$KSS.obj$dat.dim[1],ncol=x$KSS.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(x$KSS.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$KSS.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(x$KSS.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
  }else{
    if(x$KSS.obj$additive.effects=="time"){
      par(mfrow=c(1,1))
      plot.ts(x$KSS.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      par(mfrow=c(1,1))
    }
    if(x$KSS.obj$additive.effects=="twoways"){
      par(mfrow=c(1,2))
      plot.ts(x$KSS.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(x$KSS.obj$Add.Ind.Eff,each=x$KSS.obj$dat.dim[1]),
                     nrow=x$KSS.obj$dat.dim[1],ncol=x$KSS.obj$dat.dim[2]),
              main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      par(mfrow=c(1,1))
    }
    if(x$KSS.obj$additive.effects=="individual"){
      par(mfrow=c(1,1))
      matplot(matrix(rep(x$KSS.obj$Add.Ind.Eff,each=x$KSS.obj$dat.dim[1]),
                     nrow=x$KSS.obj$dat.dim[1],ncol=x$KSS.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
    par(mfrow=c(1,1))
    }
  }
}

