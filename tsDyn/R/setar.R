
#' @export
## Copyright (C) 2005,2006,2009  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY;without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

#SETAR model contructor	(sequential conditional LS)
#	th: threshold. If not specified, a grid of reasonable values is tried
#	m: general autoregressive order (mL=mH)
#	mL: autoregressive order below the threshold ('Low')
#	mH: autoregressive order above the threshold ('High')
#	nested: is this a nested call? (useful for correcting final model df)
#	trace: should infos be printed?
setar <- function(x, m, d=1, steps=d, series, mL,mM,mH, thDelay=0, mTh, thVar, th, trace=FALSE, nested=FALSE,include = c("const", "trend","none", "both"), common=c("none", "include","lags", "both"), model=c("TAR", "MTAR"), ML=seq_len(mL),MM=seq_len(mM), MH=seq_len(mH), nthresh=1,trim=0.15, type=c("level", "diff", "ADF"), restriction=c("none","OuterSymAll","OuterSymTh") ){
# 1: preliminaries
# 2:  Build the regressors matrix and Y vector
# 3: Set-up of transition variable
# 4: Search of the treshold if th not specified by user
# 5: Build the threshold dummies and then the matrix of regressors
# 6: compute the model, extract and name the vec of coeff
# 7: return the infos



###SETAR 1: preliminaries
  include<-match.arg(include)
  type<-match.arg(type)
  model<-match.arg(model)
  common<-match.arg(common)
  restriction<-match.arg(restriction)
  if(missing(m))
    m <- max(ML, MH, ifelse(nthresh==2, max(MM),0),thDelay+1)
  if(!missing(th)){
    if(length(th)==2)
      nthresh<-2
  }
  if(missing(series))
    series <- deparse(substitute(x))
  
  if(common%in%c("both", "lags")&type!="ADF")
    stop("Arg common ==", common, " only for ADF specification\n")   
  if(restriction=="OuterSymTh") 
    stop("Currently not implemented")
  if(restriction=="OuterSymAll"&nthresh==2){
    warning("With restriction ='OuterSymAll', you can only have one th. Changed to nthresh=1\n")
    nthresh<-1}
### SETAR 2:  Build the regressors matrix and Y vector
  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
  
  if(type=="level"){
    xx <- getXX(str)
    yy <- getYY(str)
  }
  else{ 
    if(type=="diff"){
      xx <- getdXX(str)
      yy <- getdYY(str)
    }
    else if(type=="ADF"){
      xx <- cbind(getdX1(str),getdXX(str))
      yy <- getdYY(str)
    }
    str$xx<-xx
    str$yy<-yy
  }
  
  externThVar <- FALSE #mmh should maybe be removed
##Lags selection
  ###ML
  if(missing(ML)) {		#ML: different lags
    if (missing(mL)) {	#mL: suit of lags
      if(missing(m))
	cat("arg m is missing")
      else
	mL <- m
      if (trace) 
	cat("Using maximum autoregressive order for low regime: mL =", m,"\n")
    }
    ML <- seq_len(mL)
  }
  
  ###MH
  if(missing(MH)) {
    if (missing(mH)) {
      mH <- m
      if (trace) 
        cat("Using maximum autoregressive order for high regime: mH =", m,"\n")
    }
    MH <- seq_len(mH)
  }
  
  ###MM
  if(missing(MM)) {
    if (missing(mM)) {
      mM <- m
      if (trace&nthresh==2) 
        cat("Using maximum autoregressive order for middle regime: mM =", m,"\n")
    }
    MM <- seq_len(mM)
  }


### SETAR 3: Set-up of transition variable (different from selectSETAR)
#two models: TAR or MTAR (z is differenced)
#three possibilitiees for thVar:
#mTh: combination of lags. Z is matrix nrow(xx) x 1
#thVar: external variable, if thDelay specified, lags will be taken, Z is matrix/vector nrow(xx) x thDelay
#former args not specified: lags of explained variable (SETAR), Z is matrix nrow(xx) x (thDelay)
# Default: thDelay=0

  maxTh<-max(thDelay)
  SeqmaxTh<-seq_len(maxTh+1)
  ###combination of lagged values
  if(!missing(mTh)) {
    if(length(mTh) != m) 
      stop("length of 'mTh' should be equal to 'm'")
    z <- xx %*% mTh #threshold variable
    dim(z) <- NULL
  } 
  ###external variable
  else if(!missing(thVar)) {
    thvarLength<-length(thVar)
    thVarTheoLen<-if(!missing(thDelay)) nrow(xx)+maxTh else nrow(xx)
    TD<-if(!missing(thDelay)) maxTh+1 else 1
    if(thvarLength>thVarTheoLen) {
      thVar <- thVar[1:thVarTheoLen]
      if(trace) 
	cat("Using only first", thVarTheoLen, "elements of thVar\n")
    }
    else if(thvarLength<thVarTheoLen) {
      if(thvarLength!=nrow(xx))
	stop("thVar has not enough/too much observations when taking thDelay")
      else{
	TD<-1
	thDelay<-0
	if(trace)
	  cat("thdelay set to default value 0")
      }
    }
    z <- embed(thVar, TD)
    externThVar <- TRUE
  }
  
  ### lagged values
  else {
    if(max(thDelay)>=m) 
      stop(paste("thDelay too high: should be < m (=",m,")"))
    if(model=="TAR"){
      if(type =="level")
	z <- getXX(str)[,SeqmaxTh]
      else
	z<-getXX(str)[-1,SeqmaxTh]
      #z2<-embedd(x, lags=c((0:(m-1))*(-d), steps) )[,1:m,drop=FALSE] equivalent if d=steps=1
      #z4<-embed(x,m+1)[,-1]
    }
    else{
      if(max(thDelay)==m-1){
	if(type =="level"){
	  z<-getdXX(str)[, SeqmaxTh]
	  xx<-xx[-1,, drop=FALSE]
	  yy<-yy[-1]
	  str$xx<-xx
	  str$yy<-yy
	}
	else 
	  z<-getdXX(str)[, SeqmaxTh]
      }
      else{
	if(type =="level")
	  z<-getdXX(str,same.dim=TRUE)[,SeqmaxTh]
	else
	  z<-getdXX(str)[,SeqmaxTh]
      }
    }
  } 

z<-as.matrix(z)

###includes const, trend (identical to selectSETAR)
  if(include=="none" && any(c(ML,MM,MH)==0))
    stop("you cannot have a regime without constant and lagged variable")
  constMatrix<-buildConstants(include=include, n=nrow(xx)) #stored in miscSETAR.R
  incNames<-constMatrix$incNames #vector of names
  const<-constMatrix$const #matrix of none, const, trend, both
  ninc<-constMatrix$ninc #number of terms (0,1, or 2)

### SETAR 4: Search of the treshold if th not specified by user
#if nthresh==1, try over a reasonable grid (30), if nthresh==2, whole values
#call the function selectSETAR
  if (missing(th)|length(thDelay)>1) { 
    search<-selectSETAR(x, m, d=d, th=list(ngrid="All"), steps=d, series, mL=mL, mH=mH,mM=mM, thDelay=thDelay, mTh, thVar, trace=trace, include = include, common=common, model=model, ML=ML,MH=MH, MM=MM,nthresh=nthresh,trim=trim,criterion = "SSR", plot=FALSE,max.iter=2, type=type, restriction=restriction)
    thDelay<-search$bests[1]
    th<-search$bests[2:(nthresh+1)]
    nested<-TRUE

    if(trace) {
      cat("Selected threshold: ", th,"\n")
      cat("Selected delay: ", thDelay,"\n")
    }
    # missing(th)<-FALSE
  }
  
  z<-z[,thDelay+1, drop=FALSE]
  
  if(restriction=="none")
    transV<-z
  else if(restriction=="OuterSymAll")
    transV<-abs(z)
  
### SETAR 5: Build the threshold dummies and then the matrix of regressors

	#check number of observations)
  if(restriction=="none"){
    if(nthresh==1){
      isL <- ifelse(z<=th, 1, 0)
      isM<-NA
      isH <- 1-isL}	
    else{
      isL <- ifelse(z<=th[1], 1, 0)
      isH <- ifelse(z>th[2], 1, 0)
      isM <- 1-isL-isH
    }
  }
  else if(restriction=="OuterSymAll"){
      isL <- ifelse(z<= -abs(th), 1, 0)
      isH <- ifelse(z>abs(th), 1, 0)
      isM <- 1-isL-isH
   }
  
  reg<-if(nthresh==1&restriction=="none") isL+2*isH else isL+2*isM+3*isH
  regime <- c(rep(NA,length(x)-length(reg)), reg)
    
  nobs<-na.omit(c(mean(isL),mean(isM),mean(isH)))	#N of obs in each regime
  if(min(nobs)<trim-0.01){
    warning("\nWith the threshold you gave (", th, ") there is a regime with less than trim=",100*trim,"% observations (",paste(round(100*nobs,2), "%, ", sep=""), ")\n", call.=FALSE)
  }
  if(min(nobs)==0)
    stop("With the threshold you gave, there is a regime with no observations!", call=FALSE)
		
	#build the X matrix
  if(type=="ADF"){
    exML<-c(1, ML+1)
    exMH<-c(1, MH+1)
    exMM<-if(nthresh==2) c(1, MM+1) else NULL
  }
  else{
    exML<-ML
    exMH<-MH
    exMM<-if(nthresh==2) MM else NULL
  }

	 
  if(nthresh==1){
    funBuild1<-switch(common, "include"=buildXth1Common, "none"=buildXth1NoCommon, "both"=buildXth1LagsIncCommon, "lags"=buildXth1LagsCommon)
    xxLH<-funBuild1(gam1=th, thDelay=0, xx=xx,trans=transV, ML=exML, MH=exMH,const, trim=trim)
  } else{
    funBuild2<-switch(common, "include"=buildXth2Common, "none"=buildXth2NoCommon, "both"=buildXth2LagsIncCommon, "lags"=buildXth2LagsCommon)
    xxLH<-funBuild2(gam1=th[1],gam2=th[2],thDelay=0,xx=xx,trans=transV, ML=exML, MH=exMH, MM=exMM,const,trim=trim)
  }

### SETAR 6: compute the model, extract and name the vec of coeff
  res <- lm.fit(xxLH, yy)
  if(any(is.na(res$coefficients)))
    warning("Problem with the regression, it may arrive if there is only one unique value in the middle regime")
#Coefficients and names
  res$coefficients <- c(res$coefficients, th)

  t<-type #just shorter
  if(common=="none"){
    if(nthresh==1){
      co<-c(getIncNames(incNames,ML), getArNames(ML,t), getIncNames(incNames,MH), getArNames(MH,t),"th")
    }else{
      co<-c(getIncNames(incNames,ML), getArNames(ML,t), getIncNames(incNames,MM), getArNames(MM,t), getIncNames(incNames,MH), getArNames(MH,t),"th1","th2")
    }
  } else if(common=="include"){
    if(nthresh==1){
      co<-c(incNames, getArNames(ML,t), getArNames(MH,t),"th")
    } else {
      co<-c(incNames, getArNames(ML,t), getArNames(MM,t), getArNames(MH,t),"th1","th2")
    }
  } else if(common=="both"){
    if(nthresh==1){
      co<-c(incNames, "phiL.1", "phiH.1", paste("Dphi.", ML),"th")
    } else {
      co<-c(incNames, "phiL.1", "phiM.1","phiH.1", paste("Dphi.", ML),"th1","th2")
    }
  } else if(common=="lags"){ #const*isL,xx[,1]*isL,xx[,1]*(1-isL),const*isH, xx[,-1]
    if(nthresh==1){
      co<-c(getIncNames(incNames,ML), "phiL.1", "phiH.1", getIncNames(incNames,ML), paste("Dphi.", ML),"th")
    } else {
      co<-c(getIncNames(incNames,ML), "phiL.1", getIncNames(incNames,MM), "phiM.1", getIncNames(incNames,MM), "phiH.1",paste("Dphi.", ML),"th")
    }
  }
  
  names(res$coefficients) <- na.omit(co)
  colnames(xxLH) <- head(na.omit(co), -nthresh)

###check for unit roots
  if(type=="level"){
    isRootH<-isRoot(res$coefficients, regime="H", lags=MH)
    isRootL<-isRoot(res$coefficients, regime="L", lags=ML)
    if(nthresh==2)
      isRootM<-isRoot(res$coefficients, regime="M", lags=MM)
  }

### SETAR 7: return the infos
  res$k <- if(nested) (res$rank+nthresh) else res$rank	#If nested, 1/2 more fitted parameter: th
  res$thDelay<-thDelay
  res$fixedTh <- if(nested) FALSE else TRUE
  res$mL <- max(ML)
  res$mH <- max(MH)
  res$mM <- if(nthresh==2) max(MH) else NULL
  res$ML <- ML
  res$MH <- MH
  res$MM <- if(nthresh==2)  MM else NULL
  res$externThVar <- externThVar
  res$thVar <- z
  res$incNames<-incNames
  res$common<-common  	#wheter arg common was given by user
  res$nthresh<-nthresh 	#n of threshold
  res$model<-model
  res$type<-type
  res$restriction<-restriction
  res$trim<-trim
  res$regime<-regime
  res$RegProp <- c(mean(isL),mean(isH))
  res$timeAttributes <- attributes(x)
	
  if(nthresh==2|restriction=="OuterSymAll")
    res$RegProp <- c(mean(isL),mean(isM),mean(isH))
  res$VAR<-as.numeric(crossprod(na.omit(res$residuals))/(nrow(xxLH)))*solve(crossprod(xxLH))
  if(!externThVar) {
    if(missing(mTh)) {
      mTh <- rep(0,m)
      mTh[thDelay+1] <- 1
    }
    res$mTh <- mTh
  }
  return(extend(nlar(str,	
	  coefficients=res$coef,
	  fitted.values=res$fitted.values,
	  residuals=res$residuals,
	  k=res$k,
	  model=data.frame(yy,xxLH),
	  model.specific=res), "setar"))
	#}
}

getSetarXRegimeCoefs <- function(x, regime=c("L","M","H")) {
  reg <- match.arg(regime)
  
  #generate vector of "^phiL|^const.L|^trend.L"
  nam1 <- paste(c("^phi", "^const.", "^trend."), reg, sep="")
  nam2 <- paste(nam1, collapse="|")
  # select that vector into coefficients
  x$coef[grep(nam2, names(x$coef))]
}



#gets a vector with names of the arg inc
getIncNames<-function(inc,ML){
  ninc<-length(inc)
  letter<-deparse(substitute(ML))
  letter<-sub("M","",letter)
  paste(inc, rep(letter,ninc), sep=".")
}

#get a vector with names of the coefficients
getArNames<-function(ML, type=c("level", "diff", "ADF")){
  phi<-ifelse(type=="level", "phi", "Dphi")
  if(any(ML==0)){
    return(NA)
  } else{
    letter<-deparse(substitute(ML))
    letter<-sub("M","",letter)
    dX1<- if(type=="ADF") paste("phi", letter, ".1", sep="") else NULL
    c(dX1, paste(paste(phi,letter, sep=""),".", ML, sep=""))
  }
}

#' @S3method print setar
print.setar <- function(x, ...) {
  NextMethod(...)
  x.old <- x
  x <- x$model.specific
  order.L <- x$mL
  order2.L <- length(x$ML)
  order.H <- x$mH
  order2.H <- length(x$MH)
  common <- x$common
  nthresh<-x$nthresh
  externThVar <- x$externThVar
  cat("\nSETAR model (",nthresh+1,"regimes)\n")
  cat("Coefficients:\n")
  if(common=="none"){
    lowCoef <- getSetarXRegimeCoefs(x.old, "L")
    highCoef<- getSetarXRegimeCoefs(x.old, "H")
    cat("Low regime:\n")
    print(lowCoef, ...)
    if(nthresh==2){
      midCoef<- getSetarXRegimeCoefs(x.old, "M")
      cat("\nMid regime:\n")
      print(midCoef, ...)}
    cat("\nHigh regime:\n")
    print(highCoef, ...)
  } else {
    print(x.old$coeff[-length(x.old$coeff)], ...)
  }
  thCoef<-getTh(coef(x.old))
  cat("\nThreshold:")
  if(x$model=="MTAR"){
    cat("\nMomentum Threshold (MTAR) Adjustment")
    D<-"D"
  }  else{
    D<-NULL
  }
  cat("\n-Variable: ")
  if(externThVar){
    cat("external")
  } else {
    cat('Z(t) = ')
    cat('+ (',format(x$mTh[1], digits=2), paste(")",D," X(t)", sep=""), sep="")
    if(length(x$mTh)>1)
      for(j in 1:(length(x$mTh) - nthresh)) {
        cat('+ (', format(x$mTh[j+1], digits=2), paste(")",D,"X(t-", j, ")", sep=""), sep="")
      }
    cat('\n')
  }
  cat("-Value:", format(thCoef, digits=4))
  if(x$fixedTh) cat(" (fixed)")
  cat("\n")
  cat("Proportion of points in ")
  if(nthresh==1&x$restriction=="none")
    cat(paste(c("low regime:","\t High regime:"), percent(x$RegProp, digits=4,by100=TRUE)), "\n")
  else
    cat(paste(c("low regime:","\t Middle regime:","\t High regime:"), percent(x$RegProp, digits=4,by100=TRUE)), "\n")
  invisible(x)
}

#' @S3method summary setar
summary.setar <- function(object, ...) {
  ans <- list()
  mod <- object$model.specific
  order.L <- mod$mL
  order2.L <- length(mod$ML)
  order.H <- mod$mH
  order2.H <- length(mod$MH)
  nthresh<-mod$nthresh		#number of thresholds
  common<-mod$common
  ans$lowCoef <- getSetarXRegimeCoefs(object, "L")
  ans$highCoef<- getSetarXRegimeCoefs(object, "H")
  ans$thCoef <- getTh(coef(object))
  ans$fixedTh <- mod$fixedTh
  ans$externThVar <- mod$externThVar
  ans$lowRegProp <- mod$lowRegProp
  n <- getNUsed(object$str)
  coef <- coef(object, hyperCoef=FALSE) #all coeffients except of the threshold
  p <- length(coef)			#Number of slope coefficients
  resvar <- mse(object) * n / (n-p)
  Qr <- mod$qr
  p1 <- 1:p
  est <- coef[Qr$pivot[p1]]
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) #compute (X'X)^(-1) from the (R part) of the QR decomposition of X.
  se <- sqrt(diag(R) * resvar) #standard errors
  tval <- est/se			# t values
  coef <- cbind(est, se, tval, 2*pt(abs(tval), n-p, lower.tail = FALSE))
  dimnames(coef) <- list(names(est), c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
  ans$coef <- coef
  ans$mTh <- mod$mTh
  extend(summary.nlar(object), "summary.setar", listV=ans)
}

#' @S3method print summary.setar
print.summary.setar <- function(x, digits=max(3, getOption("digits") - 2),
	signif.stars = getOption("show.signif.stars"), ...) {
	NextMethod(digits=digits, signif.stars=signif.stars, ...)
	cat("\nCoefficient(s):\n\n")
	printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)		
	cat("\nThreshold")
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
	cat("\nValue:", format(x$thCoef, digits=4))
	if(x$fixedTh) cat(" (fixed)")
	cat('\n')
	invisible(x)
}

#' @S3method coef setar
#Coef() method: hyperCoef=FALSE won't show the threshold coef
coef.setar <- function(object, hyperCoef=TRUE, ...){
  co <- object$coefficients
  if(!hyperCoef) co <- head(co, -length(getTh(object)))
  co
}

#' @S3method vcov setar
vcov.setar <- function(object,hyperCoef=TRUE, ...){
  mod <- object$model.specific
  
  nthresh<-mod$nthresh
  n <- length(object$str$x)
  coef <- coef(object, hyperCoef=FALSE) #all coeffients except of the threshold
  p <- length(coef)    	#Number of slope coefficients
  
  # residual variance
  resvar <- mse(object) * n / (n-p)
  
  # (X'X)^(-1)
  Qr <- mod$qr
  p1 <- 1:p
  est <- coef[Qr$pivot[p1]]
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) #compute (X'X)^(-1) from the (R part) of the QR decomposition of X.
  
  ## result
  res <- R*resvar
  
  if(hyperCoef){
    res <- cbind(rbind(res,0),0)
    if(mod$nthresh==2) res <- cbind(rbind(res,0),0)
  }
  
  return(res)
}

#' @S3method plot setar
plot.setar <- function(x, ask=interactive(), legend=FALSE, regSwStart, regSwStop, ...) {
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  par(ask=ask)
  NextMethod(ask=ask, ...)
  str <- x$str
  xx <- getXX(str)
  yy <- getYY(str)
  nms <- colnames(xx)
  m <- str$m
  d <- str$d
  lags <- c((0:(m-1))*(-d), str$steps)
  xxyy <- getXXYY(str)
  restriction<-x$model.specific$restriction
  nthresh<-x$model.specific$nthresh
  regime<-x$model.specific$regime
  x.old <- x
  x <- c(x, x$model.specific)
  series <- str$x
  z <- x$thVar
  th <- getTh(coef(x))
	  
  regime.id<-regime
  if(length(regime)<=300) {
    pch <- c(20,23)[regime.id]
    cex <- 1
  } else {
    pch <- '.'
    cex <- 4
  }
  ##Phase plot
  for(j in 1:m) {
    plot(xxyy[,j], xxyy[,m+1], xlab=paste("lag", -lags[j]), ylab=paste("lag", -lags[m+1]),
         col=regime.id, pch=pch, cex=cex, ...)
    lines.default(xxyy[,j], x.old$fitted, lty=2)
    if(nthresh==2)
      title("Curently not implemented for nthresh=2!")
    if(legend)
      legend("topleft", legend=c("low","high"), pch=pch[c(1,1)], col=1:2, merge=FALSE, title="regime")
  }
  ##Regime switching plot
  sta <- 1
  sto <- length(regime.id)
  if(!missing(regSwStart))
    sta <- regSwStart
  if(!missing(regSwStop))
    sto <- regSwStop
  t <- sta:sto
  regime.id <- regime.id[t]
  series <- series[t+(m*d)]
  ylim <- range(series)
  l <- ylim[1] * 0.9
  h <- ylim[2] * 1.1
  ylim[1] <- ylim[1] * 0.8
  ylim[2] <- ylim[2] * 1.2
  x0 <- t
  x1 <- t+1
  y0 <- series[t]
  y1 <- series[t+1]
  par(mar=c(0,4,4,0))	#number of lines of margin to be specified on the 4 sides of the plot
  plot(t, series, type=ifelse(length(regime)<=300, "p", "n"), ax=FALSE, ylab="time series values", main="Regime switching plot")
  legRSP<-c("low",if(nthresh==2|restriction=="OuterSymAll") "middle","high")
  legend("topright", legend=legRSP, pch=pch[c(1,1)], col=seq_len(length(legRSP)), merge=FALSE, title="regime", lty=1)
  abline(h=th)
  if(restriction=="OuterSymAll")
    abline(h=-th)
  axis(2)		#adds an axis on the left
  segments(x0,y0,x1,y1,col=regime.id)	#adds segments between the points with color depending on regime
  ##plot for the transition variable
  par(mar=c(5,4,4,2))
  layout(matrix(1:2, ncol=1))
  if(restriction=="OuterSymAll"){
    nthresh<-2
    th<-c(-th, th)
  }
  plot1(th, nthresh,usedThVar=x$model.specific$thVar) #shows transition variable, stored in TVARestim.R
  plot2(th, nthresh,usedThVar=x$model.specific$thVar, trim=x$model.specific$trim) #ordered transition variabl
  invisible(x)
}


oneStep.setar <- function(object, newdata, itime, thVar, ...){
  mL <- object$model.specific$mL
  mH <- object$model.specific$mH
  phi1 <- object$coefficients[1:(mL+1)]
  phi2 <- object$coefficients[mL+1+ 1:(mH+1)]
  th <- object$coefficients[mL+mH+3]
  ext <- object$model.specific$externThVar
  if(ext)	{
    z <- thVar[itime]
  } else {
    z <- newdata %*% object$model.specific$mTh
    dim(z) <- NULL
  }
  z <- (z<=th)+0
  if(nrow(newdata)>1) {
    xL <- cbind(1,newdata[,1:mL])
    xH <- cbind(1,newdata[,1:mH])
  } else {
    xL <- c(1,newdata[,1:mL])
    xH <- c(1,newdata[,1:mH])
  }
  (xL %*% phi1) * z + (xH %*% phi2) * (1-z)
}



#' Latex representation of fitted setar models
#'
#' Produce LaTeX output of the SETAR model. 
#'
#'
#' @param object fitted setar model (using \code{\link{nlar}})
#' @param digits options to be passed to \code{\link{format}} for formatting
#' numbers
#' @param label LaTeX label passed to the equation
#' @param \dots Not used
#' @author Antonio, Fabio Di Narzo
#' @seealso \code{\link{setar}}, \code{\link{nlar-methods}}
#' @keywords ts
#' @examples
#'
#' mod.setar <- setar(log10(lynx), m=2, thDelay=1, th=3.25)
#' toLatex(mod.setar)
#'
#' @method toLatex setar 
#' @S3method toLatex setar 
toLatex.setar <- function(object, digits=3, label, ...) {
  obj <- object
  mod<-obj$model.specific
  if(mod$nthresh > 1 & mod$common!="none")
    stop("currently not implemented")
  model<-mod$model
  type<-mod$type
  res <- character()
  ML <- if(type=="ADF") c(mod$ML, max(mod$ML)+1) else mod$ML #workaround...
  MH <- if(type=="ADF") c(mod$MH, max(mod$MH)+1) else mod$MH #workaround...
  steps <- obj$str$steps
  d <- obj$str$d
  betaL <- formatSignedNum(getSetarXRegimeCoefs(obj, "L"), digits=digits, ...)
  betaH <- formatSignedNum(getSetarXRegimeCoefs(obj, "H"), digits=digits, ...)
  th <- formatSignedNum(getTh(coefficients(obj)), digits=digits, ...)
  res[1] <- "\\["
  Xt<-if(type=="level") " X_{t-" else " \\Delta X_{t-"
  lab <- if(missing(label)) NULL else label
  res[2] <- paste(Xt,steps,"} = \\left\\{\\begin{array}{lr}",lab,sep="")
  translateCoefs <- function(coefs, lags, type=c("level", "diff", "ADF")) {
    ans <- ""
    if(length(lags) == (length(coefs) - 1)) { #there is a constant term
      iconst <- grep("const", names(coefs))
      stopifnot(length(iconst) == 1)
      ans <- paste(ans, coefs[iconst], " ", sep="")
      coefs <- coefs[-iconst]
    }
    m<-length(coefs)
    a<-switch(match.arg(type), "level"=m, "diff"=0, "ADF"=1)
    b<-switch(match.arg(type), "level"=0, "diff"=m, "ADF"=m-1)
    Xtt<-c(rep(" X_{t-",a), rep(" \\Delta X_{t-", b))
    for(j in seq_along(coefs)) {
      lag <- (lags[j] - 1) * d
      ans <- paste(ans, coefs[j], Xtt[j], lag, "}", sep="")
    }
    return(ans)
  }
  res[3] <- translateCoefs(betaL, ML, type=type)
  res[3] <- paste(res[3], "& Z_t \\leq ", th, "\\\\", sep="")
  res[4] <- translateCoefs(betaH, MH, type=type)
  res[4] <- paste(res[4], "& Z_t > ", th, "\\\\", sep="")
  res[5] <- "\\end{array}\\right."
  res[6] <- "\\]"
  res[7] <- ""

  if(!mod$externThVar) {
    mTh <- formatSignedNum(mod$mTh)
    m <- obj$str$m
    res[8] <- "\\["
    res[9] <- "Z_t = "
    thXt<-if(model=="TAR") " X_{t-" else " \\Delta X_{t-"
    for(j in seq_len(m)) {
      if(mod$mTh[j]==1){
        res[9] <- paste(res[9],thXt,(j-1)*d,"} ",sep="")
      }
      else if(mod$mTh[j]!=0)
        res[9] <- paste(res[9], mod$mTh[j],thXt,(j-1)*d,"} ",sep="")
    }
    res[10] <- "\\]"
    res[11] <- ""
  }
  return(structure(res, class="Latex"))
}

