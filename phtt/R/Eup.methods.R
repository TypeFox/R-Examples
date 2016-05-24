

## Methods ========================================================================================

Eup <- function(formula,
    additive.effects = c("none", "individual", "time", "twoways"),
    dim.criterion    = c("PC1", "PC2", "PC3", "BIC3", "IC1", "IC2" , "IC3", "IPC1", "IPC2", "IPC3"),
    d.max            = NULL,
    sig2.hat         = NULL,
    factor.dim       = NULL,
    double.iteration = TRUE,
    start.beta       = NULL,
    max.iteration    = 500,
    convergence      = 1e-6,
    restrict.mode    = c("restrict.factors", "restrict.loadings"),
    ...){
  UseMethod("Eup")
}

print.Eup <- function(x,...){
  cat("Call:\n")
  print(x$call)
  
  cat("\nCoeff(s) of the Observed Regressor(s) :\n\n")
  slope.para <- x$slope.para
  if(x$is.intercept){
    inter <- matrix(x$Intercept, 1, 1)
    colnames(inter) <- ""
    rownames(inter) <- "(Intercept)"
    slope.para <- rbind(signif(inter,digits=3), signif(slope.para,digits=3))
    slope.para <- signif(slope.para, 3)
  }
  print(t(slope.para))
  cat("\nAdditive Effects Type: ", as.name(x$additive.effects)," \n")
  cat("\nDimension of the Unobserved Factors:", x$used.dim," \n")
  cat("\nNumber of iterations:", x$Nbr.iteration,"\n\n")
#  cat("\nNOTE: If the panel data dimensions,'N' and  'T', are proprtional, i.e., \n") 
#  cat("the ratio 'N/T' converges asymptotically to a constant, and correlation\n")
#  cat("and/or heteroscedasticity affect the idiosyncratic errors, the user can\n") 
#  cat("specify the type of the  error structure  in the  'summary()' method to\n") 
#  cat("calculate the bias corrected estimator and obtain appropriate inference.\n")
}



coef.Eup <- function(object,...){
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

	if(object$restrict.mode[1] == "restrict.factors") exam.factors <- Ind.loadings.param
	else exam.factors <- Time.varying.ind.eff

    Var.shares.of.loadings.param      <- numeric(Factor.Dim)
    if(Factor.Dim > 0){
    Total.var.loadings.param          <- sum(apply(exam.factors,2,var))
    for(i in 1:Factor.Dim){
      Var.shares.of.loadings.param[i] <- round(var(c(exam.factors[,i]))/Total.var.loadings.param,
                    digits=4)*100
    }
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

residuals.Eup <- resid.Eup <- function(object,...){
  Residual.mat <- object$residuals
  return(Residual.mat)
}

summary.Eup <- function(object, error.type =c(1, 2, 3, 4, 5, 6, 7, 8), kernel.weights = NULL,...){
	if(!all((error.type)%in%c(1, 2, 3, 4, 5, 6, 7, 8))) 
		stop("The user specified argument 'error.type' should be on of the following numbers: 1, 2, 3, 4, 5, 6, 7, or 8.")


  ## Coefficients:
  #if(missing(kernel.weights)) kernel.weights = NULL
  TAB.obj  <-  Eup.inference(Eup.Obj=object, error.type = error.type[1], kernel.weights = kernel.weights)
  TAB  <- TAB.obj$inf.result
  yy <- sum(diag(crossprod(object$orig.Y - mean(object$orig.Y))))

  if(error.type[1]==7|error.type[1]==8){
	BC.pca <- pca.fit(dat = TAB.obj$BC.remeind, given.d = object$used.dim, restrict.mode = object$restrict.mode)
	Resids <- TAB.obj$BC.remeind - BC.pca$fitted.values
	nr <- ncol(Resids)
	sig2.hat <- sum(diag(var(Resids)))*(nr-1)/object$degrees.of.freedom
	}
  else {
	BC.pca <- NULL
  	Resids <- object$residuals
	sig2.hat <- object$sig2.hat
	}

  Res.outpt <- signif((summary(as.vector(Resids))), digits=3)[-4]

  names(Res.outpt) <- c("Min", "1Q", "Median", "3Q", "Max")
  ee <- sum(c(Resids)^2)
  R2 <- signif(1 - ee/yy, 4)
  
 
  ## Result:
  result        <- list(Res.outpt    = Res.outpt,
                        coefficients = signif(TAB, 3),
                        R2 		 = R2,
				sig2.hat     = sig2.hat,
                        Eup.obj      = object, 
				error.type   = error.type[1], BC.pca = BC.pca)                
  class(result) <- "summary.Eup"
  result
}

print.summary.Eup <- function(x,...){

  ## Call
  cat("Call:\n")
  print(x$Eup.obj$call)
  ## Residuals:
  cat("\nResiduals:\n")
  print(x$Res.outpt)
  cat("\n")
  ## Beta-Coeffs
  if(x$error.type == 7 | x$error.type == 8) cat("\n Bias Corrected Slope-Coefficients:\n")
  else cat("\n Slope-Coefficients:\n")
  printCoefmat(x$coefficients)
  cat("\nAdditive Effects Type: ", as.name(x$Eup.obj$additive.effects)," \n")
  cat("\nDimension of the Unobserved Factors:", x$Eup.obj$used.dim," \n")
  cat("\nResidual standard error:", signif(sqrt(x$sig2.hat), 4), "on",
            x$Eup.obj$degrees.of.freedom, "degrees of freedom, ", "\nR-squared:", x$R2,"\n")
}


coef.summary.Eup <- function(object,...){
    if(object$Eup.obj$is.intercept){
    Intercept <- object$coefficients[1,1]
    Slope.Coef <- as.matrix(object$coefficients[-1,1])
    }
    else{
    Intercept <- 0
    Slope.Coef <- as.matrix(object$coefficients[,1])
    }
    
    
    if(object$Eup.obj$additive.effects== "individual"| object$Eup.obj$additive.effects== "twoways")
    Add.Ind.Eff <- object$Eup.obj$ColMean[,1, drop = FALSE] - object$Eup.obj$ColMean[,-1, drop = FALSE]%*%Slope.Coef - c(Intercept)
    else Add.Ind.Eff <- rep.int(0, object$Eup.obj$dat.dim[2])

    if(object$Eup.obj$additive.effects== "time"| object$Eup.obj$additive.effects== "twoways")
    Add.Tim.Eff <- object$Eup.obj$RowMean[,1, drop = FALSE] - object$Eup.obj$RowMean[,-1, drop = FALSE]%*%Slope.Coef - c(Intercept)
    else Add.Tim.Eff <- rep.int(0, object$Eup.obj$dat.dim[1])

    if(is.null(object$BC.pca)){
    Common.factors <- object$Eup.obj$unob.factors  
    Ind.loadings.param <- object$Eup.obj$ind.loadings   
    Time.varying.ind.eff <- object$Eup.obj$unob.fact.stru
    }
    else{
    Common.factors <- object$BC.pca$factors
    Ind.loadings.param <- object$BC.pca$loadings
    Time.varying.ind.eff <- object$BC.pca$fitted.values
    }
    Factor.Dim <- object$Eup.obj$used.dim

	if(object$Eup.obj$restrict.mode[1] == "restrict.factors") exam.factors <- Ind.loadings.param
	else exam.factors <- Time.varying.ind.eff

    Var.shares.of.loadings.param      <- numeric(Factor.Dim)
    if(Factor.Dim > 0){
    Total.var.loadings.param          <- sum(apply(exam.factors,2,var))
    for(i in 1:Factor.Dim){
      Var.shares.of.loadings.param[i] <- round(var(c(exam.factors[,i]))/Total.var.loadings.param,
                    digits=4)*100
    }
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

residuals.summary.Eup <- resid.Eup <- function(object,...){
  if(is.null(object$BC.pca))
  Residual.mat <- object$Eup.obj$residuals
  else Residual.mat <- object$BC.pca$orig.values - object$BC.pca$fitted.values
  return(Residual.mat)
}


plot.summary.Eup <- function(x,...){
  if(coef(x)$Factor.Dim == 0 & x$Eup.obj$additive.effects=="none"){
    stop("Neither an estimated factor structure nor additive effects to plot.")
  }
  if(coef(x)$Factor.Dim > 0){
    if(x$Eup.obj$additive.effects=="none"){
      par(mfrow=c(1,2))
      matplot(coef(x)$Common.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(coef(x)$Time.varying.ind.eff,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,3))
      plot.ts(coef(x)$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(coef(x)$Common.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(coef(x)$Time.varying.ind.eff,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,4))
      plot.ts(coef(x)$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(coef(x)$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(coef(x)$Common.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(coef(x)$Time.varying.ind.eff,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,3))
      matplot(matrix(rep(coef(x)$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(coef(x)$Common.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="o",...)
      matplot(coef(x)$Time.varying.ind.eff,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
  }else{
    if(x$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,1))
      plot.ts(coef(x)$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,2))
      plot.ts(coef(x)$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(coef(x)$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
              main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,1))
      matplot(matrix(rep(coef(x)$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
    par(mfrow=c(1,1))
    }
  }
}

