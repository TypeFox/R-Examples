gaus <-
function(data=data, B=B, taus, formula, basis=NULL, alpha=0.05, var, load, rearrange=F, rearrange.vars="quantile", 
                 uniform=F, se="unconditional", average=T, nderivs=1, method="fn")
{
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)] 
  mf$drop.unused.levels <- TRUE 
  mf[[1]] <- as.name("model.frame") 
  mf <- eval.parent(mf) 
  mt<-attr(mf,"terms")
  y<-model.response(mf,"numeric")
  x<-model.matrix(mt,mf)
  if(all(x[,1]==rep(1,dim(x)[1]))){x <- x[,-1]} # Remove intercept from x matrix
  
  
  z<-as.vector(data[,var]) # Save vector with variable of interest
  f<-formula(y~x) # Save variables locally and create formula
  length.taus=length(taus)
  basis.type<-"noBasis"
  
  if (is.basis(basis)){
    if (basis$type=="bspline" | basis$type=="fourier"){
      basis.type <- "fda"
    } else {
      stop("basis must be a bspline basis, fourier basis, polynomial basis, or factor variable")
    }
  } else if (is.factor(basis)) {
    basis.type <- "factor"
  } else if (attr(basis,"class")[1]=="poly") {
    basis.type <- "polynomial"
  } else {
    stop("basis must be a bspline basis, fourier basis, polynomial basis, or factor variable")
  }
  
  if (rearrange==T & (nderivs!=0 | average!=F)) {
    stop( "Rearrangement may only be used when nderivs=0 and average=F") }
  
  if (rearrange==T & rearrange.vars!="quantile" & rearrange.vars!="var" & rearrange.vars!="both") {
    stop( "rearrange.vars must be quantile, var, or both")}
  
  if (average==T & nderivs==0) {
    stop( "This method is not available for average=T and nderivs=0")}
  
  if (average!=T & se=="unconditional") {
    stop( "Unconditional standard error estimation is only valid when average=T")
  }
  if (is.null(basis)) {
    stop( "Please input a basis")
  }
  if(basis.type=="factor" & (average!=F | nderivs!=0)){
    stop( "Fully saturated indicator series approximation is only available when average=F and nderivs=0")
  }
  if((basis.type=="polynomial") & nderivs>2) {
    stop("Estimation using polynomials is not available for nderivs>2")
  }
  
  if (basis.type=="fda") {
    basvecs<-eval.basis(z,basis)[,-1] # Create vectors of nonparametric bases evaluated at observed variable values
  } else {
    basvecs<-basis
  }
  
  f<-update.formula(f, ~ . + basvecs) # Include series regressors in formula; supress intercept (included in nonparametric basis)
  nObs <- dim(x)[1]
  
  # Generate loading vector if not input by user
  if (is.null(load)) {
    
    modelMatrix <- model.matrix(f, data)
    nVars <- dim(modelMatrix)[2]-dim(basvecs)[2]-1 # number of model variables, excluding the intercept and the var of interest
    
    
    if (average==T) { # generate loading vector
      
      # Get derivative of formula, not including the nonparametric components
      parDeriv <- formulaDeriv(formula,var,data=data,nDerivs=nderivs)   
      if (basis.type=="fda") {
        load <- cbind(matrix(apply(parDeriv, 2, mean), nrow=1), matrix(apply(as.matrix(eval.basis(z,basis,Lfdobj=nderivs)[,-1]), 2, mean), nrow=1)) 
      } else {
        load <- cbind(matrix(apply(parDeriv, 2, mean), nrow=1), matrix(apply(as.matrix(poly.wrap(x=z,degree=max(attr(basis,"degree")),nderivs=nderivs,coefs=attributes(basis)$coefs)), 2, mean), nrow=1)) 
      }
    } else {
      
      if (nderivs==0) {
        load <- as.matrix(aggregate(modelMatrix, by=list(z), FUN=load.sum)[,-1])
      } else {
        # Get derivative of formula, not including the nonparametric components
        parDeriv <- formulaDeriv(formula,var,data=data,nDerivs=nderivs)
        if (basis.type=="fda"){
          load <- as.matrix(aggregate(cbind(parDeriv, eval.basis(z,basis,Lfdobj=nderivs)[,-1]), by=list(z), FUN=mean)[,-1])
        } else {
          load <- as.matrix(aggregate(cbind(parDeriv, poly.wrap(x=z,degree=max(attr(basis,"degree")),nderivs=nderivs,coefs=attributes(basis)$coefs)), by=list(z), FUN=mean)[,-1])
        }
      }  
    }
    
  }
  
  Q <- msqrt(t(modelMatrix) %*% modelMatrix)
  z.unique <- as.matrix(aggregate(z, by=list(z), FUN=mean)[,-1]) # generate vector with unique values of z
  loadLength <- dim(load)[1]
  loadWidth <- dim(load)[2]
  if (basis.type=="fda") {
    nVars <- dim(modelMatrix)[2]-length(basis$names)
  } else if (basis.type=="polynomial"){
    nVars <- dim(modelMatrix)[2]-(max(attr(basis,"degree"))+1)
  } else {
    nVars <- dim(modelMatrix)[2]-nlevels(basis)
  }
  length.taus=length(taus)
  
  if (average==T) {   
    if (basis.type=="fda") {
      covMat <- cov(cbind(matrix(0, nrow=nObs, ncol=nVars+1), eval.basis(z,basis,Lfdobj=nderivs)[,-1])) #to be used in unconditional standard errors
    } else {
      covMat <- cov(cbind(matrix(0, nrow=nObs, ncol=nVars+1), poly.wrap(x=z,degree=max(attr(basis,"degree")),nderivs=nderivs,coefs=attributes(basis)$coefs))) #to be used in unconditional standard errors
    }
  } else {     
    covMat <- NULL       
  }
  
  gaussian	<- list(qfits=NULL, point.est=NULL, var.unique=NULL, CI=NULL, CI.oneSided=NULL, std.error=NULL, pvalues=NULL, load=NULL) # set up output vector
  gaussian$qfits <- vector("list",length.taus)   
  
  Gstat <- NULL
  
  A  	<- kronecker(as.matrix(taus), matrix(1, ncol = length.taus ));
  C		<- kronecker(matrix(1,nrow=length.taus), t(taus));
  BB.sigma	<- (A <= C) * A + (A > C) * C - A*C; # generate sigma for brownian bridge
  
  BB  <- t(rmnorm(n = B, mean = rep(0, length.taus*ncol(Q)), varcov = kronecker(BB.sigma, diag(1,nrow = ncol(Q), ncol = ncol(Q))) ) );
        
  Gstat <- array(0,dim=c(length.taus,B,loadLength))
  Hinvs <- vector("list",length.taus)
  
  err<-F
  for (i in 1:length.taus) {
    gaussian$qfits[[i]]<-tryCatch(expr=rq(f,tau=taus[i],data,method=method),error=function(e) {
      err<-T
      fit<-rq(f,tau=taus[i],data)
      return(fit)
    })
    Hinvs[[i]] <- summary.rq(gaussian$qfits[[i]], se= "ker", covariance = T)$Hinv
  }
  
  if (err==T) {
    message("Warning: rq method chosen produced an error; using method=br when necessary")
  }
  
  for (j in 1:loadLength) {
     for (i in 1:length.taus) {      
        Gstat[i, ,j]  <- load[j,] %*% Hinvs[[i]] %*% Q %*% BB[((i-1)*ncol(Q) + 1):(i*ncol(Q)), ];                          
     }               
  }

  # Generate output, including t statistics
  fittedValues <- matrix(0, nrow=loadLength, ncol=length.taus) # predicted values from quantile regressions
  
  for (i in 1:length.taus){ # predict fitted values, evaluated at loading vector
    fittedValues[,i] <- load %*% gaussian$qfits[[i]]$coef
  }
  if (rearrange==T) { # perform rearrangement if requested
    if (rearrange.vars=="both") {
      if (length(taus)==1){
        fittedValues <- matrix(rearrangement(list(z.unique),as.vector(fittedValues)),nrow=loadLength,ncol=1)
      } else {
        fittedValues <- rearrangement(list(z.unique,taus),fittedValues)
      }
    } else if (rearrange.vars=="quantile") {
      if (length(taus)==1){
        fittedValues <- fittedValues
      } else {
        for (j in 1:loadLength) {
          fittedValues[j,] <- rearrangement(list(taus),as.vector(fittedValues[j,]))
        }
      }
    } else if (rearrange.vars=="var") {
      for (i in 1:length.taus) {
        fittedValues[,i] <- rearrangement(list(z.unique),as.vector(fittedValues[,i]))
      }
    }
  }
  
  Covs         <- array (0, dim=c(length.taus, loadWidth,loadWidth)) # covariance matrix for each tau
  Tn           <- array (0, dim=c(loadLength,length.taus,B)) # array of t stats
  Tn.oneSided  <- array (0, dim=c(loadLength,length.taus,B, 2)) # array of t stats for one-sided CIs
  Kn.oneSided  <- array (0, dim=c(loadLength,length.taus,2)) # vector of alpha and 1-alpha sample quantile of one sided t stats
  CI           <- array (0, dim=c(loadLength,length.taus,2)) # two vectors of CIs at each tau (two-sided CI)
  CI.oneSided  <- array (0, dim=c(loadLength,length.taus,2)) # two vectors of CIs at each tau (one-sided CIs)
  Kn           <- matrix(0, nrow=loadLength, ncol=length.taus) # matrix of 1-alpha sample quantile of t stat
  stdError     <- matrix(0, nrow=loadLength, ncol=length.taus) # estimate of standard error for each tau at each data point
  tDenom       <- matrix(0, nrow=loadLength, ncol=length.taus) # denominator of t stats
  
  for (i in 1:length.taus){ # generate estimates of standard errors

    Covs[i, ,] <- summary.rq(gaussian$qfits[[i]], se = "ker", covariance = T)$cov
    for (j in 1:loadLength){
      if (se=="unconditional") {
        stdError[j,i] <- sqrt((t(load[j,]) %*% Covs[i, ,] %*% load[j,]) + t(gaussian$qfits[[i]]$coef) %*% covMat %*% gaussian$qfits[[i]]$coef / nObs)
      } else { 
        stdError[j,i] <- sqrt(t(load[j,]) %*% Covs[i, ,] %*% load[j,])
      }
    }
      
    for (j in 1:loadLength){ # generate t statistics
      tDenom[j,i] <- stdError[j,i]
      Tn[j,i, ] <- abs(Gstat[i, ,j] / tDenom[j,i])
      Tn.oneSided[j,i, ,1] <- Gstat[i, ,j] / tDenom[j,i]
      Tn.oneSided[j,i, ,2] <- Gstat[i, ,j] / tDenom[j,i]
    }
  }
  
  if (uniform==T) { # if uniform inference requested, make t statistics uniform
    for (b in 1:B){
      Tn[ , ,b] <- matrix(max(Tn[ , ,b]),nrow=loadLength,ncol=length.taus)
      Tn.oneSided[, ,b,1] <- matrix(max(Tn.oneSided[ , ,b,1]),nrow=loadLength,ncol=length.taus)
      Tn.oneSided[, ,b,2] <- matrix(max(Tn.oneSided[ , ,b,2]),nrow=loadLength,ncol=length.taus)  
    }
  }
  
  for (j in 1:loadLength){ # generate K (quantiles of t stats)
    for (i in 1:length.taus) {
      Kn[j,i] <- quantile(Tn[j,i, ],probs=(1-alpha),na.rm=T)
      Kn.oneSided[j,i,1] <- quantile(Tn.oneSided[j,i, ,1], probs=(1-alpha),na.rm=T)
      Kn.oneSided[j,i,2] <- quantile(Tn.oneSided[j,i, ,2], probs=alpha,na.rm=T)
      
      # generate confidence intervals
      CI[j,i,1] <- fittedValues[j,i] - Kn[j,i]*stdError[j,i]
      CI[j,i,2] <- fittedValues[j,i] + Kn[j,i]*stdError[j,i]
      CI.oneSided[j,i,1] <- fittedValues[j,i] - Kn.oneSided[j,i,1]*stdError[j,i]
      CI.oneSided[j,i,2] <- fittedValues[j,i] + Kn.oneSided[j,i,1]*stdError[j,i]
    } 
  }
  
  k.onesided <- max(fittedValues/stdError)
  k.star.onesided <- apply((Gstat-aperm(array(fittedValues,dim=c(loadLength,length.taus,B)),perm=c(2,3,1)))/aperm(array(stdError,dim=c(loadLength,length.taus,B)),perm=c(2,3,1)), c(1,2), max, na.rm=T)
  pValue.onesided.neg <- mean(k.star.onesided >= k.onesided)
  pValue.onesided.pos <- mean(k.star.onesided <= k.onesided)
  
  
  k.twosided <- max(abs(fittedValues/stdError))
  k.star.twosided <- apply(abs((Gstat-aperm(array(fittedValues,dim=c(loadLength,length.taus,B)),perm=c(2,3,1)))/aperm(array(stdError,dim=c(loadLength,length.taus,B)),perm=c(2,3,1))), c(1,2), max, na.rm=T)
  pValue.twosided <- mean(k.star.twosided >= k.twosided)
  
  pvalues <- c(pValue.onesided.neg, pValue.onesided.pos, pValue.twosided)

  gaussian$var.unique <- z.unique
  
  gaussian$point.est <- matrix(NA, nrow=loadLength, ncol=length.taus) 
  gaussian$point.est <- fittedValues
  
  gaussian$CI <- CI
  gaussian$CI.oneSided <- CI.oneSided
  gaussian$std.error <- stdError
  gaussian$pvalues <- pvalues
  gaussian$load <- load
  
  return(gaussian);
}
