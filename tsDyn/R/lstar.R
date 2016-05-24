
#' @export
## Copyright (C) 2005, 2006, 2007/2006  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

# LSTAR fitter
#	mTh: as in setar
#	phi1, phi2, c, gamma: initial guesses for model parameters
#	trace: should infos be printed?
#	control: 'control' options to be passed to optim
lstar <- function(x, m, d=1, steps=d, series, mL, mH, mTh, thDelay,
                  thVar, th, gamma, trace=TRUE, include = c("const", "trend","none", "both"), control=list(), starting.control=list())
{

  include<-match.arg(include)

  if(missing(m))
    m <- max(mL, mH, thDelay+1)

  if(missing(series))
    series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)

  xx <- str$xx
  yy <- str$yy
  
  externThVar <- FALSE
  if (missing(mL)) {
    mL <- m
    if (trace) 
      cat("Using maximum autoregressive order for low regime: mL =", m,"\n")
  }
  if (missing(mH)) {
    mH <- m
    if (trace) 
      cat("Using maximum autoregressive order for high regime: mH =", m,"\n")
  }
  if(!missing(thDelay)) {
    if(thDelay>=m) 
      stop(paste("thDelay too high: should be < m (=",m,")"))
    z <- xx[,thDelay+1]
  }
  else if(!missing(mTh)) {
    if(length(mTh) != m) 
      stop("length of 'mTh' should be equal to 'm'")
    z <- xx %*% mTh #threshold variable
    dim(z) <- NULL
  }
  else if(!missing(thVar)) {
    if(length(thVar) > nrow(xx)) {
      z <- thVar[1:nrow(xx)]
      if(trace) 
        cat("Using only first", nrow(xx), "elements of thVar\n")
    }
    else {
      z <- thVar
    }
    externThVar <- TRUE
  }
  else {
    if(trace) 
      cat("Using default threshold variable: thDelay=0\n")
    z <- xx[,1]
    thDelay = 0
  }
  ## Build regressors matrix
  constMatrix<-buildConstants(include=include, n=nrow(xx)) #stored in miscSETAR.R
  incNames<-constMatrix$incNames #vector of names
  const<-constMatrix$const #matrix of none, const, trend, both
  ninc<-constMatrix$ninc #number of terms (0,1, or 2)

  xxL <- cbind(const,xx[,1:mL])
  xxH <- cbind(const,xx[,1:mH])

  #Fitted values, given parameters
  #phi1: vector of 'low regime' parameters
  #phi2: vector of 'high regime' parameters
  #g: smoothing parameter
  #c: threshold value
  #Model covariates are 'xxL', 'xxH' and 'x', as defined in the
  #   beginning of that function
  F <- function(phi1, phi2, g, th, type=1){
    if(type==1){
      xxL %*% phi1 + (xxH %*% phi2) * G(z, g, th)
    } else {
      (xxL %*% phi1)* (1-G(z, g, th)) + (xxH %*% phi2) * G(z, g, th)
    }
  }

  F_bind <- function(xxL, xxH, g, th, type=1){
    if(type==1){
      cbind(xxL , xxH * G(z, g, th))
    } else {
      cbind(xxL * (1- G(z, g, th)), xxH * G(z, g, th))
    }
    
  }
#Automatic starting values####################
  if(missing(th) || missing(gamma)) {
    if (trace)
      cat("Performing grid search for starting values...\n");

    bestCost <- Inf;

    # Set list of defaults:
    start.con<-list(
		    nTh=200, 
		    trim=0.1,
		    nGamma=40,
		    gammaInt=c(1,100), 
		    thInt=NA, 
		    candidates=NA
    )
    # Add if user defined, check if names correspond (code taken from optim)
    nmsC <- names(start.con)
    start.con[(namc <- names(starting.control))] <- starting.control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in starting.control: ", paste(noNms, collapse = ", "))

  ## Set grid search values
    interv.Th <- quantile(as.ts(z), c(start.con$trim, 1-start.con$trim)) # "trim" percentil of z
    if(any(is.na(start.con$thInt))) interv.Th <- c(min(start.con$thInt[1], interv.Th[1], na.rm=TRUE), min(start.con$thInt[2], interv.Th[2], na.rm=TRUE))
    Gammas <- seq(start.con$gammaInt[1], start.con$gammaInt[2], length.out=start.con$nGamma)
    ths <- seq(interv.Th[1], interv.Th[2], length.out=start.con$nTh) 

    IDS <- as.matrix(expand.grid(Gammas, ths) )

    if(!is.na(start.con$candidates)){
      li <- start.con$candidates
      if(length(li)!=2 | any(names(li)!=c("th", "gamma")) | length(li[[1]])!=length(li[[2]])){
	stop("Error in specification of starting.control$candidates: should be a list with element 'th' and 'gamma' of same length\n")
      }
      IDS <- rbind(IDS, cbind(li[["gamma"]], li[["th"]]))
    }
  ## Grid search: Loop over values
    for(i in 1:nrow(IDS)){

      # We fix the linear parameters.
      cost <- crossprod(lm.fit(F_bind(xxL, xxH, g=IDS[i,1], th=IDS[i,2]), yy)$residuals)

      if(cost <= bestCost) {
	bestCost <- cost;
	gamma <- IDS[i,1]
	th <- IDS[i,2]
      }
    }

    if (trace) {
      cat("Starting values fixed: gamma = ", gamma,", th = ", th, 
          "; SSE = ", bestCost, "\n")
      if(gamma%in%start.con$gammaInt) cat("Grid search selected lower/upper bound gamma (was: ", start.con$gammaInt, "]). 
					  Might try to widen bound with arg: 'starting.control=list(gammaInt=c(1,200))'\n")
    }
  }

  # Fix the linear parameters one more time
#   new_phi<- lm.fit(cbind(xxL, xxH * transFun(z, gamma, th)), yy)$coefficients
#   phi1 <- new_phi[1:(mL+1)]
#   phi2 <- new_phi[(mL+2):(mL + mH + 2)]

  # Computes the gradient 
  #
  # Returns the gradient with respect to the error
  gradEhat <- function(p)
    {
      gamma <- p[1]  #Extract parms from vector p
      th    <- p[2] 	     #Extract parms from vector p
      new_phi<- lm.fit(cbind(xxL, xxH * G(z, gamma, th)), yy)$coefficients
      phi1 <- new_phi[1:(mL+ninc)]
      phi2 <- new_phi[(mL+ninc+1):(mL + mH + 2*ninc)]

      y.hat <- F(phi1, phi2, gamma, th)
      e.hat <- yy - y.hat

      fX <- sigmoid(gamma * (z - th));
      dfX <- dsigmoid(fX);
      
      gGamma <- as.vector(xxH %*% phi2) * as.vector(dfX * (z - th));
      gTh <-        - as.vector(xxH %*% phi2) * as.vector(gamma * dfX);

      J = - cbind(gGamma, gTh) / sqrt(str$n.used)
      
      return(2 * t(e.hat) %*% J)
      
    }
  
  #Sum of squares function
  #p: vector of parameters
  SS <- function(p) {
    gamma <- p[1]   #Extract parms from vector p
    th <- p[2]      #Extract parms from vector p

    trans <- G(z, gamma, th)
    m_trans <- mean(trans, na.rm=TRUE)
    pen <- if(min(m_trans, 1-m_trans, na.rm=TRUE)< 0.05) 1/(0.05-m_trans) else 0

    # First fix the linear parameters
    xx <- F_bind(xxL, xxH, g=gamma, th=th)
    if(any(is.na(as.vector(xx)))) {
      message('lstar: missing value during computations')
      return (Inf)
    }
    crossprod(lm.fit(xx, yy)$residuals) + pen
  }
 
  ## Numerical minimization##########
  p <- c(gamma, th)   #pack parameters in one vector
  res <- optim(p, SS, gradEhat, hessian = FALSE, method="BFGS", control = control)
  if(trace){
    if(res$convergence!=0){
      if(res$convergence==1) {
	cat("Convergence problem code 1. You might want to increase maximum number of iterations by setting 'control=list(maxit=1000)'\n")
      } else {
	cat("Convergence problem. Convergence code: ",res$convergence,"\n")
      }
    } else {
      cat("Optimization algorithm converged\n")
    }
  }
  phi_2<- lm.fit(F_bind(xxL, xxH, g=res$par[1],th= res$par[2]), yy)$coefficients
  coefnames_L <- c(if(ninc>0) paste(incNames,"L",sep="."), paste("phiL", 1:mL, sep="."))
  coefnames_H <- c(if(ninc>0) paste(incNames,"H",sep="."), paste("phiH", 1:mH, sep="."))
  coefnames <- c(coefnames_L, coefnames_H)

  names(phi_2) <-coefnames
  names(res$par) <- c("gamma", "th")

  ## Optimization: second quick step to get hessian for all parameters ########
  SS_2 <- function(p) {
    phi1 <- p[coefnames_L]	#Extract parms from vector p
    phi2 <- p[coefnames_H]	#Extract parms from vector p
    y.hat <- F(phi1, phi2, g=p["gamma"], th=p["th"])
    crossprod(yy - y.hat)
  }
  res$par <- c(phi_2,res$par)
  if(as.numeric(R.Version()$minor)<15){
    res$hessian <- optim(res$par, SS_2, method="L-BFGS-B", hessian=TRUE)$hessian
  } else {
    res$hessian <- optimHess(res$par , SS_2)
  }

  if(trace & qr(res$hessian, 1e-07)$rank != length(res$par)){
    cat("Problem: singular hessian\n")
  }


  # Results storing ################
  coefs <- res$par
  names(coefs) <- c(coefnames,"gamma", "th") 
  gamma <- coefs["gamma"]
  th  <- coefs["th"]
  if (trace) cat("Optimized values fixed for regime 2 ",
                 ": gamma = ", gamma, ", th = ", th,"; SSE = ", res$value, "\n");
  
  res$coefficients <- coefs
  res$mL <- mL
  res$mH <- mH
  res$externThVar <- externThVar
  if(!externThVar) {
    if(missing(mTh)) {
      mTh <- rep(0,m)
      mTh[thDelay+1] <- 1
    }
    res$mTh <- mTh
  }
  res$thVar <- c(rep(NA, length(x)-length(z)),z)
  res$fitted <- F(coefs[coefnames_L], 
                  coefs[coefnames_H], gamma, th, type=1)
  res$residuals <- yy - res$fitted
  dim(res$residuals) <- NULL	#this should be a vector, not a matrix
  res$k <- length(res$coefficients)
  res$ninc<-ninc
  res$include<-include
  res$timeAttributes <- attributes(x)

  mod.return <- data.frame(yy,F_bind(xxL, xxH, g=gamma, th=th))
  colnames(mod.return) <- c("yy", coefnames)

################################

  return(extend(nlar(str, 
                     coefficients=res$coef,
                     fitted.values =res$fitted,
                     residuals =res$residuals,
                     k   =res$k,
		     model = mod.return,
                     model.specific=res),
                "lstar"))
}

	

#############################################
  #Transition function G: moved to star.R

#' @S3method print lstar
print.lstar <- function(x, ...) {
  NextMethod(...)
  cat("\nLSTAR model\n")
  x <- x$model.specific
  ninc<-x$ninc
  order.L <- x$mL
  order.H <- x$mH
  lowCoef <- x$coef[grep("phiL|const\\.L|trend\\.L", names(x$coef))]
  highCoef <- x$coef[grep("phiH|const\\.H|trend\\.H", names(x$coef))]
  gammaCoef <- x$coef["gamma"]
  thCoef <- x$coef["th"]
  externThVar <- x$externThVar
  
  cat("Coefficients:\n")
  cat("Low regime:\n")
  print(lowCoef, ...)
  cat("\nHigh regime:\n")
  print(highCoef, ...)
  cat("\nSmoothing parameter: gamma =", format(gammaCoef, digits=4),"\n")
  cat("\nThreshold")
  cat("\nVariable: ")
  if(externThVar)
    cat("external")
  else {
    cat('Z(t) = ')
    cat('+ (',format(x$mTh[1], digits=2), ') X(t) ', sep="")
    if(length(x$mTh)>1)
      for(j in 1:(length(x$mTh) - 1)) {
        cat('+ (', format(x$mTh[j+1], digits=2), ') X(t-', j, ')', sep="")
      }
    cat('\n')
  }
  cat("\nValue:", format(thCoef, digits=4), "\n")
  invisible(x)
}

#' @S3method summary lstar
summary.lstar <- function(object, ...) {
  ans <- list()  

  ## SE for ML estimates (from optim), taken from 'arma.R' in package tseries
  coef<- object$coefficients
  n <- object$str$n.used
  vc <- vcov(object)
  if(all(is.na(vc))) {
    se <- rep(NA, length(coef))
  } else {
    se <- sqrt(diag(vc))
  }

  tval <- coef/ se
  coefMat<- cbind(coef, se, tval, 2 * ( 1 - pnorm( abs(tval) ) ) )
  dimnames(coefMat) <- list(names(coef), c(" Estimate"," Std. Error"," t value","Pr(>|z|)"))
  ans$coefficients<-coefMat

  ## Non-linearity test############
  xx <- object$str$xx
  sX <- tail(object$model.specific$thVar, -(object$str$n.used-nrow(object$str$xx)))
  dim(sX) <- NULL
  xx1<- xx*sX		#predictors set B (approximated non-linear component)
  yy <- object$str$yy
  modA <- lm(yy ~ xx)
  modB <- update(modA, . ~ . + xx1)
  a <- anova(modA, modB)
  ans$nlTest.value <- a[["F"]][2]
  ans$nlTest.pval  <- a[["Pr(>F)"]][2]
  
###############################
  ninc    <- object$model.specific$ninc
  order.L <- object$model.specific$mL
  order.H <- object$model.specific$mH
  ans$lowCoef <- object$coef[1:(order.L+ninc)]
  ans$highCoef<- object$coef[(order.L+1+ninc):(order.H+2*ninc)]
  ans$thCoef <- object$coef["th"]
  ans$externThVar <- object$model.specific$externThVar
  ans$mTh <- object$model.specific$mTh
  return(extend(summary.nlar(object), "summary.lstar", listV=ans))
}

#' @S3method print summary.lstar
print.summary.lstar <- function(x, digits=max(3, getOption("digits") - 2),
                       signif.stars = getOption("show.signif.stars"), ...)
{
  NextMethod(digits=digits, signif.stars=signif.stars, ...)
  cat("\nCoefficient(s):\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  cat("\nNon-linearity test of full-order LSTAR model against full-order AR model\n")
  cat(" F =", format(x$nlTest.value, digits=digits),"; p-value =", format(x$nlTest.pval, digits=digits),"\n")
  cat("\nThreshold ")
  cat("\nVariable: ")
  if(x$externThVar)
    cat("external")
  else {
    cat('Z(t) = ')
    cat('+ (',format(x$mTh[1], digits=2), ') X(t) ', sep="")
    if(length(x$mTh)>1)
      for(j in 1:(length(x$mTh) - 1)) {
        cat('+ (', format(x$mTh[j+1], digits=2), ') X(t-', j, ')', sep="")
      }
    cat('\n')
  }
  invisible(x)
}

#' @S3method plot lstar
plot.lstar <- function(x, ask=interactive(), legend=FALSE,
                       regSwStart, regSwStop, ...) {
  
  op <- par(no.readonly=TRUE)
  par(ask=ask)
  NextMethod(ask=ask, ...)
  str <- x$str
  xx <- str$xx
  yy <- str$yy
  nms <- colnames(xx)
  z <- tail(x$model.specific$thVar, -(x$str$n.used-nrow(x$str$xx)))
  z <- plogis(z, x$coefficients["th"], 1/x$coefficients["gamma"])
  regime.id <- cut(z, breaks=quantile(z, 0:5/5), include.lowest=TRUE)
  regime.id <- as.numeric(regime.id)
  if(length(regime.id)<=300) {
    pch <- regime.id
    cex <- 1
  }
  else {
    pch <- '.'
    cex <- 4
  }
  palette <- rgb(0:4/4,0,0)
  for(j in 1:x$str$m) {
    plot(xx[,j], yy, xlab=nms[j], ylab=paste("lag",x$str$steps),
         col=palette[regime.id], pch=pch, cex=cex, ...)
    lines.default(xx[,j], x$model.specific$fitted, lty=2)
    if(legend) {
      labels <- c("[0;0.2]","(0.2,0.4]","(0.4;0.6]","(0.6;0.8]","(0.8;1]")
      legend("topleft", legend=labels, pch=sort(unique(regime.id)),
             col=palette[sort(unique(regime.id))], title="regime quantiles")
    }
  }
  sta <- 1
  sto <- length(regime.id)
  if(!missing(regSwStart))
    sta <- regSwStart
  if(!missing(regSwStop))
    sto <- regSwStop
  t <- sta:sto
  regime.id <- regime.id[t]
  m <- x$str$m
  d <- x$str$d
  series <- x$str$x[t+(m*d)]
  ylim <- range(series)
  l <- ylim[1] * 0.9
  h <- ylim[2] * 1.1
  ylim[1] <- ylim[1] * 0.8
  ylim[2] <- ylim[2] * 1.2
  x0 <- t
  x1 <- t+1
  y0 <- series[t]
  y1 <- series[t+1]
  par(mar=c(0,4,4,0))
  plot(t, series, type="n", ax=FALSE, ylab="time series values",
       main="Regime switching plot")
  axis(2)
  segments(x0,y0,x1,y1,col=palette[regime.id])
  par(op)
  invisible(x)
}

#' @S3method coef lstar
#Coef() method: hyperCoef=FALSE won't show the threshold/slope coef
coef.lstar <- function(object, hyperCoef=TRUE, ...){
  co <- object$coefficients
  if(!hyperCoef) co <- head(co, -2)
  co
}

#' @S3method vcov lstar
vcov.lstar <- function(object, ...){
  n <- object$str$n.used
  coef<- object$coefficients
  rank <- qr(object$model.specific$hessian, 1e-07)$rank

  if(rank != length(coef(object))) {
    vc <- matrix(NA, length(coef),length(coef))
    warning("singular Hessian\n")
  } else{
    vc <- 2*object$model.specific$value/n*solve(object$model.specific$hessian)
    if(any(diag(vc) < 0))
      warning("Hessian negative-semi definite\n")
  }

return(vc)
}

#' @S3method confint lstar
confint.lstar <- function(object, parm, level = 0.95, ...){
  confint.default(object, parm=parm, level=level, ...)
}

oneStep.lstar <- function(object, newdata, itime, thVar, ...){
  include <- object$model.specific$include
  if(!include %in%c("none","const")) stop("oneStep currently only implemented for include==const/none\n")
  ninc<-object$model.specific$ninc
  mL <- object$model.specific$mL
  mH <- object$model.specific$mH
  coefs<-object$coefficients
  phi1 <- coefs[grep("const\\.L|trend\\.L|phiL", names(coefs))]
  phi2 <- coefs[grep("const\\.H|trend\\.H|phiH", names(coefs))]
  gamma <- coefs["gamma"]
  c <- coefs["th"]
  ext <- object$model.specific$externThVar

  if(ext) {
    z <- thVar[itime]
  }
  else {
    z <- newdata %*% object$model.specific$mTh
    dim(z) <- NULL
  }
  z <- plogis(z, c, 1/gamma)

  if(nrow(newdata)>1) {
    reg<- switch(include,"const"=1, "none"=NULL)
    xL <- cbind(reg,newdata[,1:mL])
    xH <- cbind(reg,newdata[,1:mH])
  } else {
    reg<- switch(include,"const"=1, "none"=NULL)
    xL <- c(reg,newdata[,1:mL])
    xH <- c(reg,newdata[,1:mH])
  }

  xL %*% phi1 + (xH %*% phi2) * z
}


