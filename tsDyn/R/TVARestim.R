#'Multivariate Threshold Autoregressive model
#'
#'Estimate a multivariate Threshold VAR
#'
#'For fixed \code{th} and threshold variable, the model is linear, so
#'estimation can be done directly by CLS (Conditional Least Squares). The
#'search of the parameters values is made upon a grid of potential values. So
#'it is pretty slow.
#'
#'nthresh=1: estimation of one threshold model (two regimes) upon a grid of
#'\var{ngrid} values (default to ALL) possible thresholds and delays values.
#'
#'nthresh=2: estimation of two thresholds model (three regimes) Conditional on
#'the threshold found in model where nthresh=1, the second threshold is
#'searched. When both are found, a second grid search is made with 30 values
#'around each threshold.
#'
#'nthresh=3: DOES NOT estimate a 3 thresholds model, but a 2 thresholds model
#'with a whole grid over the thresholds parameters (so is really slow) with a
#'given delay, is there rather to check the consistency of the method nthresh=2
#'
#'@aliases TVAR OlsTVAR
#'@param data time series
#'@param lag Number of lags to include in each regime
#'@param include Type of deterministic regressors to include
#'@param model Whether the transition variable is taken in levels (TAR) or
#'difference (MTAR)
#'@param commonInter Whether the deterministic regressors are regime specific
#'(commonInter=FALSE) or not.
#'@param nthresh Number of thresholds
#'@param thDelay 'time delay' for the threshold variable (as multiple of
#'embedding time delay d) PLEASE NOTE that the notation is currently different
#'to univariate models in tsDyn. The left side variable is taken at time t, and
#'not t+1 as in univariate cases.
#'@param mTh combination of variables with same lag order for the transition
#'variable. Either a single value (indicating which variable to take) or a
#'combination
#'@param thVar external transition variable
#'@param trim trimming parameter indicating the minimal percentage of
#'observations in each regime
#'@param ngrid number of elements of the grid, especially for \code{nthresh=3}
#'@param gamma prespecified threshold values
#'@param around The grid search is restricted to \var{ngrid} values around this
#'point. Especially useful for \code{nthresh=3}.
#'@param plot Whether a plot showing the results of the grid search should be
#'printed
#'@param dummyToBothRegimes Whether the dummy in the one threshold model is
#'applied to each regime or not.
#'@param trace should additional infos be printed out?
#'@param trick type of R function called: \code{for} or \code{mapply}
#'@param max.iter Number of iterations for the algorithm
#'@return An object of class TVAR, with standard methods.
#'@author Matthieu Stigler
#'@seealso \code{\link{lineVar}} for the linear VAR/VECM,
#'\code{\link{TVAR.LRtest}} to test for TVAR, \code{\link{TVAR.sim}} to
#'simulate/bootstrap a TVAR.
#'@references Lo and Zivot (2001) "Threshold Cointegration and Nonlinear
#'Adjustment to the Law of One Price," Macroeconomic Dynamics, Cambridge
#'University Press, vol. 5(4), pages 533-76, September.
#'@keywords ts
#'@export
#'@examples
#'
#'data(zeroyld)
#'
#'data<-zeroyld
#'
#'TVAR(data, lag=2, nthresh=2, thDelay=1, trim=0.1, mTh=1, plot=TRUE)
#'
#'##The one threshold (two regimes) gives a value of 10.698 for the threshold and 1 for the delay. 
#' #Conditional on this values, the search for a second threshold (three regimes) gives 8.129. 
#' #Starting from this values, a full grid search finds the same values and confims 
#' #the first step estimation. 
#'
TVAR <- function(data, lag, include = c( "const", "trend","none", "both"), model=c("TAR", "MTAR"), commonInter=FALSE, nthresh=1,thDelay=1, mTh=1,thVar, trim=0.1,ngrid, gamma=NULL,  around, plot=FALSE, dummyToBothRegimes=TRUE, trace=TRUE, trick="for", max.iter=2){
y <- as.matrix(data)
Torigin <- nrow(y) 	#Size of original sample
T <- nrow(y) 		#Size of start sample
p <- lag
t <- T-p 		#Size of end sample
k <- ncol(y) 		#Number of variables
t<-T-p			#Size of end sample


if(is.null(colnames(data)))
	colnames(data)<-paste("Var", c(1:k), sep="")
if(max(thDelay)>p)
	stop("Max of thDelay should be smaller or equal to the number of lags")
if(dummyToBothRegimes==FALSE&nthresh!= 1) 
	warning("The 'dummyToBothRegimes' argument is only relevant for one threshold models")
model<-match.arg(model)
include<-match.arg(include)

Y <- y[(p+1):T,] #
Z <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

if(include=="const")
	Z<-cbind(1, Z)
else if(include=="trend")
	Z<-cbind(seq_len(t), Z)
else if(include=="both")
	Z<-cbind(rep(1,t),seq_len(t), Z)
if(commonInter & include!="const")
	stop("commonInter argument only avalaible with include = const")
npar <- ncol(Z)			#Number of parameters


########################
### Threshold variable
########################

###External threshold variable
if (!missing(thVar)) {		
        if (length(thVar) > Torigin) {
		z <- thVar[seq_len(Torigin)]
		warning("The external threshold variable is not of same length as the original variable")
        }
        else
		z <- thVar
	z<-embed(z,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
combin<-NULL
} ###Combination (or single value indicating position) of contemporaneous variables
else {
	if (!length(mTh)%in%c(1,k))
		stop("length of 'mTh' should be equal to the number of variables, or just one")
	if(length(mTh)==1) {
		if(mTh>k)
			stop("Unable to select the ",mTh, "variable for the threshold. Please see again mTh ")
		combin <- matrix(0,ncol=1, nrow=k)
		combin[mTh,]<-1
	}
	else 
		combin<-matrix(mTh,ncol=1, nrow=k)
	zcombin <- y %*% combin
	if(model=="MTAR"){
		if(max(thDelay)<p)
			z<-embed(diff(zcombin),p)[,seq_len(max(thDelay))+1]
		else if(max(thDelay)==p){
			z<-embed(diff(zcombin),p+1)[,seq_len(max(thDelay))+1]
			z<-rbind(0,as.matrix(z))}
	}
	else
		z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
}

trans<-as.matrix(z)

###############################
###Grid for transition variable
###############################

allgammas <- sort(unique(trans[,1]))
nga <- length(allgammas)
ninter <- round(trim*nga)
gammas <- allgammas[(trim*nga):((1-trim)*nga)]



if(!missing(ngrid)){
	gammas <- allgammas[seq(from=ceiling(trim*nga), to=floor((1-trim)*nga), length.out=ngrid)]
}
if(!missing(gamma)){
	gammas<-gamma
	plot<-FALSE
}

Y_t<-t(Y)					#dim k x t

if(!missing(around)){
	if(missing(ngrid)) ngrid<-20
	if(length(around)==1)
		gammas <- aroundGrid(around, allgammas,ngrid,trim, trace=trace)
	if(length(around)==2) {
		gammas <- aroundGrid(around[1], allgammas,ngrid,trim, trace=trace)
		gammas2 <- aroundGrid(around[2], allgammas,ngrid,trim, trace=trace)
	}
}

######################
###One threshold functions					
######################

#Model with dummy applied to only one regime
loop1_onedummy <- function(gam1, thDelay){
	##Threshold dummies
	dummyDown <- ifelse(trans[,thDelay]<=gam1, 1,0) * Z
	ndown<-mean(dummyDown)
	regimeDown<-dummyDown*Z
	##SSR
	if(min(ndown, 1-ndown)>=trim){
		Z1 <- cbind(regimeDown, Z)		# dim t x k(p+1) 
		res <- crossprod(c(lm.fit(x=Z1, y=Y)$resid))
	}	else {
		res<-NA
	}
	return(res)
} #end of the function


#Model with dummy applied to both regimes
loop1_twodummy <- function(gam1, thDelay){
	##Threshold dummies
	d1<-ifelse(trans[,thDelay]<=gam1, 1,0)
	ndown<-mean(d1)
	##SSR
	if(min(ndown, 1-ndown)>=trim){
	  Z1 <- cbind(d1 * Z, (1-d1)*Z)    # dim k(p+1) x t
	  res <- crossprod(c(lm.fit(x=Z1, y=Y)$resid))
	}	else{
		res<-NA
	}
	return(res)
} #end of the function

#Model with dummy applied to both regimes and a common intercept
loop1_twodummy_oneIntercept <- function(gam1, thDelay){
	##Threshold dummies
	d1<-ifelse(trans[,thDelay]<=gam1, 1,0)
	ndown<-mean(d1)
	if(min(ndown, 1-ndown)>=trim){
		Z1 <- cbind(1,d1 * Z[,-1], (1-d1)*Z[,-1])		# dim k(p+1) x t
		res <- crossprod(c(lm.fit(x=Z1, y=Y)$resid)) 
	} else {
		res<-NA
	}
	return(res)
} #end of the function


#######################
###Two thresholds functions
#######################

loop2 <- function(gam1, gam2,thDelay){
	##Threshold dummies
	dummydown <- ifelse(trans[,thDelay]<=gam1, 1, 0)
	regimedown <- dummydown*Z
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[,thDelay]>gam2, 1, 0)
	regimeup <- dummyup*Z
	nup <- mean(dummyup)
	##SSR from TVAR(3)
	#print(c(ndown,1-nup-ndown,nup))
	if(min(nup, ndown, 1-nup-ndown)>trim){
		Z2 <- cbind(regimedown, (1-dummydown-dummyup)*Z, regimeup)		# dim k(p+1) x t
		res <- crossprod(c(lm.fit(x=Z2, y=Y)$resid)) 
	}
	else
		res <- NA
	return(res)
}

loop2_oneIntercept <- function(gam1, gam2,thDelay){
	##Threshold dummies
	dummydown <- ifelse(trans[,thDelay]<=gam1, 1, 0)
	regimedown <- dummydown*Z[,-1]
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[,thDelay]>gam2, 1, 0)
	regimeup <- dummyup*Z[,-1]
	nup <- mean(dummyup)
	##SSR from TVAR(3)
	#print(c(ndown,1-nup-ndown,nup))
	if(min(nup, ndown, 1-nup-ndown)>trim){
		Z2 <- cbind(1,regimedown, (1-dummydown-dummyup)*Z, regimeup)		# dim k(p+1) x t	
		res <- crossprod(c(lm.fit(x=Z2, y=Y)$resid)) 
	}
	else
		res <- NA
	return(res)
}
############################
###Search for one threshold
############################
if(!missing(around))
	gammas <- aroundGrid(around,allgammas,ngrid,trim, trace=trace)
if(dummyToBothRegimes){
	if(commonInter)
		func<-loop1_twodummy_oneIntercept
	else
		func <-loop1_twodummy}
else	
	func <- loop1_onedummy


bestone <- onesearch(thDelay,gammas, fun=func, trace=trace, gamma=gamma, plot=plot)
bestThresh <- bestone$bestThresh
bestDelay <- bestone$bestDelay
allThSSR<-bestone$allres


############################
###Search for two threshold
############################


if(nthresh==2){
  ###Conditionnal step
  if(commonInter)
	func2<-loop2_oneIntercept
  else
	func2<-loop2
# secondBestThresh<-condiStep(allgammas, threshRef=bestThresh, delayRef=bestDelay,ninter=ninter, fun=func2)$newThresh
# step2FirstBest<-condiStep(allgammas, threshRef=secondBestThresh, delayRef=bestDelay,ninter=ninter, fun=func2)

last<-condiStep(allgammas, threshRef=bestThresh, delayRef=bestDelay, fun=func2, trim=trim, trace=trace)
i<-1
while(i<max.iter){
	b<-condiStep(allgammas, threshRef=last$newThresh, delayRef=bestDelay, fun=func2, trim=trim, trace=trace)
	if(b$SSR<last$SSR){	#minimum still not reached
		i<-i+1
		last<-b}
	else{			#minimum reached
		i<-max.iter
		last<-b}
}

bests<-c(last$threshRef, last$newThresh)

###Alternative step: grid around the points from first step
smallThresh <- min(bests)		#bestThresh,secondBestThresh)
gammasDown <- aroundGrid(around=smallThresh,allgammas,ngrid=30, trim=trim, trace=trace)

bigThresh <- max(bests)			#bestThresh,secondBestThresh)
gammasUp <- aroundGrid(around=bigThresh,allgammas,ngrid=30, trim=trim, trace=trace)

bestThresh<-grid(gammasUp, gammasDown, fun=func2, method=trick, thDelay=bestDelay, trace=trace)

}#end if nthresh=2

###Search both thresholds with d given

if(nthresh==3){
bestDelay <- thDelay
if(missing(gamma)==FALSE){
	gammas <- gamma[1]
	gammas2 <- gamma[2]
	ninter<- 2
	cat("To be corrected!!")
}
if(missing(around)==FALSE){
	if(length(around)!=2)
		stop("Please give two thresholds possible values to search around")
	gammas <- aroundGrid(min(around), allgammas, ngrid=ngrid, trim=trim, trace=trace)
	gammas2 <- aroundGrid(max(around), allgammas, ngrid=ngrid, trim=trim, trace=trace)
}
else {
	gammas2 <- gammas
	if(length (gammas) * length(gammas2)/2>10000)
		cat("The function will compute about", length(gammas)*length(gammas2)/2, "operations. Take a coffee and come back\n")
}
if(length(thDelay)>1) stop("length of thDelay should not be bigger than 1. The whole search is made only upon the thresholds with given delay")

store3 <- matrix(NA,ncol=length(gammas2), nrow=length(gammas))

###Loop
for(i in seq_len(length(gammas))){
	gam1 <- gammas[i]
	for (j in seq_len(length(gammas))){
		gam2 <- gammas2[j]
		store3[i,j] <- loop2(gam1, gam2, thDelay=bestDelay)
	}
}

position <- which(store3==min(store3, na.rm=TRUE), arr.ind=TRUE)
r <- position[1]
c <- position[2]

gamma1 <- gammas[r]
gamma2 <- gammas2[c]
bestThresh <- c(gamma1, gamma2)

}#end n

#############
###Best Model
#############
if(commonInter)
	val<- 1
else
	val<- -(seq_len(ncol(Z)))

if(nthresh==1){
	dummydown <- ifelse(trans[,bestDelay]<=bestThresh, 1, 0)
	ndown <- mean(dummydown)
	regimeDown <- dummydown*Z[,-val]
	dummyup<-1-dummydown
	if(dummyToBothRegimes) 
		regimeUp<-dummyup*Z[,-val]
	else regimeUp<-Z
	if(commonInter)
		Zbest<-cbind(1,regimeDown,regimeUp)
	else
		Zbest <- cbind(regimeDown,regimeUp)		# dim k(p+1) x t
}

if(nthresh==2|nthresh==3){
	dummydown <- ifelse(trans[,bestDelay]<=bestThresh[1], 1,0)
	ndown <- mean(dummydown)
	regimedown <- dummydown*Z[,-val]
	dummyup <- ifelse(trans[,bestDelay]>bestThresh[2], 1,0)
	nup <- mean(dummyup)
	regimeup <- dummyup*Z[,-val]
	dummymid<-1-dummydown-dummyup
	if(commonInter)
		Zbest <- cbind(1,regimedown,dummymid*Z[,-1], regimeup)	# dim k(p+1) x t
	else
		Zbest <- cbind(regimedown,dummymid*Z, regimeup)	# dim k(p+1) x t
}

Zbest_t <- t(Zbest)
reg<-if(nthresh==1) dummydown+2*dummyup else dummydown+2*dummymid+3*dummyup
regime <- c(rep(NA, T-t), reg)


final <- lm.fit(x=Zbest, y=Y)

Bbest <- t(final$coef)
fitted <- final$fitted
resbest <- final$residuals

SSRbest <- as.numeric(crossprod(c(resbest)))
nparbest<-nrow(Bbest)*ncol(Bbest)

Sigmabest<-matrix(1/t*crossprod(resbest),ncol=k,dimnames=list(colnames(data), colnames(data)))
SigmabestOls<-Sigmabest*(t/(t-ncol(Bbest)))

###naming and dividing B
rownames(Bbest) <- paste("Equation", colnames(data))
LagNames<-c(paste(rep(colnames(data),p), -rep(seq_len(p), each=k)))

Bnames<-c(switch(include, const="Intercept", trend="Trend", both=c("Intercept","Trend"), none=NULL),LagNames)
Blist<-nameB(mat=Bbest, commonInter=commonInter, Bnames=Bnames, nthresh=nthresh, npar=npar)
BnamesVec<-if(class(Blist)=="list") c(sapply(Blist, colnames)) else colnames(Blist)
colnames(Bbest)<-BnamesVec

##number of obs in each regime
if(nthresh==1)
	nobs <- c(ndown=ndown, nup=1-ndown)
else if (nthresh==2)
	nobs <- c(ndown=ndown, nmiddle=1-nup-ndown,nup=nup)

###Y and regressors matrix
tZbest<-Zbest
naX<-rbind(matrix(NA, ncol=ncol(tZbest), nrow=p), tZbest)
YnaX<-cbind(data, naX)
BlistMod<-nameB(mat=Bbest, commonInter=commonInter, Bnames=Bnames, nthresh=nthresh, npar=npar,sameName=FALSE )
BnamesVecMod<-if(class(BlistMod)=="list") c(sapply(BlistMod, colnames)) else colnames(BlistMod)
colnames(YnaX)<-c(colnames(data),BnamesVecMod)

###elements to return
specific<-list()
specific$allgammas<-allgammas
specific$gammas<-gammas
specific$thDelay<-bestDelay
specific$Thresh<-bestThresh
specific$nthresh<-nthresh
specific$transCombin<-combin
specific$regime<-regime
specific$nreg<-nthresh+1
specific$nrowB<-npar
specific$nobs<-nobs
specific$oneMatrix<-commonInter
specific$threshEstim<-ifelse(is.null(gamma), TRUE, FALSE)
specific$allThSSR<-allThSSR#SSR values for all the th computed
specific$Bnames<-Bnames
specific$timeAttributes <- attributes(data[,1])

z<-list(coefficients=Blist, coeffmat=Bbest, residuals=resbest, model=YnaX, nobs_regimes=nobs, k=k, t=t, T=T,nparB=nparbest, fitted.values=fitted, lag=lag, include=include,model.specific=specific, usedThVar=trans[,bestDelay], trim=trim)
class(z)<-c("TVAR","nlVar")
attr(z, "levelTransVar")<-model
attr(z, "transVar")<-if(!missing(thVar)) "external" else "internal"
attr(z, "varsLevel")<-"level"

if(plot){
  layout(matrix(1:ifelse(z$model.specific$threshEstim,3,2), ncol=1))
  plot1(bestThresh, nthresh,usedThVar=z$usedThVar)
  plot2(bestThresh, nthresh,usedThVar=z$usedThVar, trim=z$trim)
  if(z$model.specific$threshEstim)
    plot3(bestThresh, nthresh,allTh=z$model.specific$allThSSR)
}
return(z)
}	#end of the whole function


#' @S3method print TVAR
print.TVAR<-function(x,...){
# 	NextMethod(...)
	cat("Model TVAR with ", x$model.specific$nthresh, " thresholds\n\n")
	print(x$coefficients)
	cat("\nThreshold value")
	print(paste(x$model.specific$Thresh, collapse=" "))
}

#' @S3method summary TVAR
summary.TVAR<-function(object,...){
	x<-object
	xspe<-x$model.specific
	k<-x$k
	t<-x$t
	p<-x$lag
	Z<-t(as.matrix(tail.matrix(x$model[,-c(1:k)],t)))
	###Stdev, VarCov
	Sigmabest<-matrix(1/x$t*crossprod(x$residuals),ncol=k)
	SigmabestOls<-Sigmabest*(x$t/(x$t-ncol(x$coeffmat)))
	VarCovB<-solve(tcrossprod(Z))%x%SigmabestOls
	StDevB<-matrix(diag(VarCovB)^0.5, nrow=k)
	Tvalue<-x$coeffmat/StDevB
	StDevB<-nameB(StDevB,commonInter=xspe$oneMatrix, Bnames=xspe$Bnames, nthresh=xspe$nthresh, npar=xspe$nrowB)
	Pval<-pt(abs(Tvalue), df=(ncol(Z)-nrow(Z)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(ncol(Z)-nrow(Z)), lower.tail=TRUE)
	Pval<-nameB(Pval,commonInter=xspe$oneMatrix, Bnames=xspe$Bnames, nthresh=xspe$nthresh, npar=xspe$nrowB)
	x$coefficients<-asListIfMat(x$coefficients)
	x$StDev<-asListIfMat(StDevB)
	x$Pvalues<-asListIfMat(Pval)
	x$Tvalues<-Tvalue
	x$VarCov<-asListIfMat(VarCovB)
	ab<-list()
	symp<-list()
	stars<-list()
	for(i in 1:length(x$Pvalues)){
		symp[[i]] <- symnum(x$Pvalues[[i]], corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," "))
		stars[[i]]<-matrix(symp[[i]], nrow=nrow(x$Pvalues[[i]]))
		ab[[i]]<-matrix(paste(x$coefficients[[i]],"(", x$StDev[[i]],")",stars[[i]], sep=""), nrow=nrow(x$StDev[[i]]))
		dimnames(ab[[i]])<-dimnames(x$coefficients[[1]])
	}
	attributes(ab)<-attributes(x$coefficients)
	x$stars<-stars
	x$starslegend<-symp[[1]]
	x$bigcoefficients<-ab
	x$aic<-AIC.nlVar(x)
	x$bic<-BIC.nlVar(x)
	x$SSR<-deviance(x)
	class(x)<-c("summary.TVAR", "TVAR", "nlVar")
	return(x)
}

#' @S3method print summary.TVAR
print.summary.TVAR<-function(x,digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"),...){
	coeftoprint<-list()
	for(i in 1:length(x$bigcoefficients)){
		a<-myformat(x$coefficients[[i]], digits)
		b<-myformat(x$StDev[[i]], digits)
		if(getOption("show.signif.stars"))
			stars<-x$stars[[i]]	
		else
			stars<-NULL
		coeftoprint[[i]]<-matrix(paste(a,"(", b,")",stars, sep=""), nrow=nrow(x$StDev[[1]]))
		dimnames(coeftoprint[[i]])<-dimnames(x$coefficients[[1]])
	}
	cat("Model TVAR with ", x$model.specific$nthresh, " thresholds\n")
	cat("\nFull sample size:",x$T, "\tEnd sample size:", x$t) 
	cat("\nNumber of variables:", x$k,"\tNumber of estimated parameters:", x$npar,"+",x$model.specific$nthresh)
	cat("\nAIC",x$aic , "\tBIC", x$bic, "\t SSR", x$SSR,"\n\n")
	print(noquote(coeftoprint))
	if (signif.stars) 
	        cat("---\nSignif. codes: ", attr(x$starslegend, "legend"), "\n")
	cat("\nThreshold value:",x$model.specific$Thresh)
	if(!x$model.specific$threshEstim)
		cat(" (user specified)")
	cat("\nPercentage of Observations in each regime:", percent(x$model.specific$nobs,3,TRUE), "\n")
}




plot1<-function(th,nthresh,usedThVar){
plot.ts(usedThVar, ylab="")
  title("Threshold variable used")
  abline(h=th, col=2:(nthresh+1))
  legend("topleft", lty=1, bg="white", col=2:(nthresh+1), legend=c(paste("th", 1:nthresh)))
}

plot2<-function(th,nthresh, usedThVar,trim){
nga <- length(usedThVar)
  orderedThVar<-sort(usedThVar)
  lt<-rep(".",nga)
  lt[(trim*nga):((1-trim)*nga)]<-1
  #plot(seq_len(nga), allgammas, pch=lt, xlab="Sorted values", ylab="")
  if(nthresh==1)
    numTh<-which.min(abs(orderedThVar-th))
  else{
    numTh1<-which.min(abs(orderedThVar-th[1]))
    numTh2<-which.min(abs(orderedThVar-th[2]))
    numTh<-c(numTh1, numTh2)}
  ts.plot(orderedThVar, xlab="", ylab="")
  title("Ordered threshold variable")
  abline(v=c(trim*nga,(1-trim)*nga), col=4)
  points(numTh, th, col=c(2:(nthresh+1)),cex=2)
  leg<-c(paste("trim=", trim),paste("th", 1:nthresh))
  legend("topleft", lty=1, col=c(4,2:(nthresh+1)), legend=leg, bg="white")
}

#Plot for the grid search
plot3<-function(th,nthresh, allTh){
    allDelay<-unique(allTh[,1])
    col <- rep(allDelay,length.out=nrow(allTh))+1
    if(nthresh==1)
      posBestTh<-which(allTh[,2]==th)
    else{
      posBestTh1<-which(allTh[,2]==th[1])
      posBestTh2<-which(allTh[,2]==th[2])
      posBestTh<-c(posBestTh1,posBestTh2)}
    plot(allTh[,2], allTh[,3], col=col,xlab="Threshold Value",ylab="SSR")
    title("Results of the grid search")
    points(th,allTh[posBestTh,3], col=c(2:(nthresh+1)), cex=2)
    leg<-c(paste("Threshold Delay", allDelay),(paste("th", 1:nthresh)))
    legend("topleft", pch=1, legend=leg, col=c(allDelay+1,c(2:(nthresh+1))), bg="white")
}

#' @S3method plot TVAR
plot.TVAR<-function(x,ask=interactive(), ...){
  th<-x$model.specific$Thresh
  nthresh<-x$model.specific$nthresh
  op <- par(no.readonly=TRUE)
  #par(ask=ask)
  layout(matrix(1:ifelse(x$model.specific$threshEstim,3,2), ncol=1))
  plot1(th, nthresh,usedThVar=x$usedThVar)
  plot2(th, nthresh,usedThVar=x$usedThVar, trim=x$trim)
  if(x$model.specific$threshEstim)
    plot3(th, nthresh,allTh=x$model.specific$allThSSR)
  par(op)
}


#' @S3method toLatex TVAR
toLatex.TVAR<-function(object,..., digits=4, parenthese=c("StDev","Pvalue")){
	x<-object
	parenthese<-match.arg(parenthese)
	x$coefficients<-asListIfMat(x$coefficients)
	if(inherits(x,"summary.TVAR")){
		coeftoprint <-list()
		for(i in 1:length(x$coefficients)){
			a<-myformat(x$coefficients[[i]],digits, toLatex=TRUE)
			if(parenthese=="StDev")
				b<-myformat(x$StDev[[i]],digits,toLatex=TRUE)
			else if(parenthese=="Pvalue")
				b<-myformat(x$Pvalues[[i]],digits,toLatex=TRUE)
			if(getOption("show.signif.stars"))
				stars<-paste("^{",x$stars[[i]],"}", sep="")
			else
				stars<-NULL
			coeftoprint[[i]]<-matrix(paste(a,"(",b,")",stars, sep=""),ncol=ncol(a), nrow=nrow(a))
		}#end for
	}#end if

	else{
		coeftoprint <-rapply(x$coefficients,myformat, digits=digits,toLatex=TRUE, how="list")}
	varNames<-rownames(coeftoprint[[1]])
	res<-character()
	res[1]<-"%insert in the preamble and uncomment the line you want for usual /medium /small matrix"
	res[2]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\begin{pmatrix}}{\\end{pmatrix}} %USUAL"
	res[3]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\left(\\begin{smallmatrix}}{\\end{smallmatrix}\\right)} %SMALL"
	res[4]<-"%\\usepackage{nccmath} \\newenvironment{smatrix}{\\left(\\begin{mmatrix}}{\\end{mmatrix}\\right)} %MEDIUM"
	res[5]<-"\\begin{equation}"
	res[6]<- "\\begin{smatrix} %explained vector"
	res[7]<-TeXVec(paste("X_{t}^{",seq(1, x$k),"}", sep=""))
	res[8]<- "\\end{smatrix}="
	#if(!x$model.specific$oneMatrix)
	res[length(res)+1]<- "\\left\\{"
 	res[length(res)+1]<-"\\begin{array}{ll}"
	Th<-x$model.specific$Thresh
	nthresh<-length(Th)
	###Condition for the threshold
	if(nthresh%in%c(1,2)){
		cond<-paste(c("& \\text{if Th}<","& \\text{if Th}>"), Th)}
	if(nthresh==2){
		cond[3]<-cond[2]
		cond[2]<-paste("& \\text{if }",Th[1], "< \\text{Th} <", Th[2])
 		}
	###Adds the const/trend and lags
	for(i in 1:(nthresh+1)){
		if(x$model.specific$oneMatrix){
			regimei<-coeftoprint[[1]]
			j<-i}
		else{
			regimei<-coeftoprint[[i]]
			j<-1}
 		res<-include(x, res, regimei)
		ninc<-switch(x$include, "const"=1, "trend"=1,"none"=0, "both"=2)
		res<-LagTeX(res, x, regimei, skip=ninc+x$lag*x$k*(j-1))
		res[length(res)+1]<- paste(cond[i], "\\\\")
	}
	res[length(res)+1]<-"\\end{array}"
	res[length(res)+1]<-"\\right."
	res[length(res)+1]<-"\\end{equation}"
	res<-gsub("slash", "\\", res, fixed=TRUE)
	res<-res[res!="blank"]
	return(structure(res, class="Latex"))

}


nameB<-function(mat,commonInter, Bnames, nthresh, npar, model=c("TVAR","TVECM"), TVECMmodel="All", sameName=TRUE){
  model<-match.arg(model)
  addRegLetter<-if(sameName) NULL else  c("L ", if(nthresh==1) NULL else "M ", "H ")
  if(model=="TVAR")
    sBnames<-Bnames[-which(Bnames=="Intercept")]
  else if(model=="TVECM")
    sBnames<-Bnames[-which(Bnames=="ECT")]

##1 threshold
  if(nthresh==1){
    if(commonInter){
      if(model=="TVAR"){
        colnames(mat)<-c("Intercept",paste(rep(addRegLetter, each=length(sBnames)),rep(sBnames,2)),sep="")
      } else if(model=="TVECM"){
        colnames(mat)<-c("ECT-","ECT+", sBnames)
      }
      Blist<-mat
    }else{
      colnames(mat) <- paste(rep(addRegLetter, each=length(Bnames)),rep(Bnames,2),sep="")
      Bdown <- mat[,c(1:npar)]
      Bup <- mat[,-c(1:npar)]
      Blist <- list(Bdown=Bdown, Bup=Bup)
    }
##2 thresholds
  } else{ 
    if(commonInter){
      if(model=="TVAR")
        colnames(mat)<-c("Intercept",paste(rep(addRegLetter, each=length(sBnames)),rep(sBnames,3)),sep="")
      else if(model=="TVECM")
        colnames(mat)<-c("ECT-","ECT+", sBnames)
      Blist<-mat}
    else{
      colnames(mat)<-paste(rep(addRegLetter, each=length(Bnames)),rep(Bnames,3),sep="")
      Bdown <- mat[,c(1:npar)]
      Bmiddle <- mat[,c(1:npar)+npar]
      Bup <- mat[,c(1:npar)+2*npar]		
      colnames(Bmiddle) <- Bnames
      Blist <- list(Bdown=Bdown, Bmiddle=Bmiddle,Bup=Bup)}
  }
  return(Blist)
}

onesearch <- function(thDelay,gammas, fun, trace, gamma, plot){
  grid1 <- expand.grid(thDelay,gammas)				#grid with delay and gammas
  store<-mapply(fun,thDelay=grid1[,1],gam1=grid1[,2])
  posBestThresh <- which(store==min(store, na.rm=TRUE), arr.ind=TRUE)[1]
  
  if(trace)
    cat("Best unique threshold", grid1[posBestThresh,2],"\n")
  if(length(thDelay)>1&trace)
    cat("Best Delay", grid1[posBestThresh,1],"\n")
  res<-list(bestThresh=grid1[posBestThresh,2],bestDelay=grid1[posBestThresh,1], allres=cbind(grid1,store))
  return(res)
}#end of function one search

condiStep<-function(allTh, threshRef, delayRef, fun, trim, trace=TRUE, More=NULL){
  allThUniq <- unique(allTh)
  ng <- length(allTh)
  down<-ceiling(trim*ng)
  
   #correction for case with few unique values
  if(allTh[down]==allTh[down+1]){
    sames<-which(allTh==allTh[down])
    down <-sames[length(sames)]+1
  }
  up<-floor(ng*(1-trim))
  #correction for case with few unique values
  if(allTh[up]==allTh[up-1]){
    up<-which(allTh==allTh[up])[1]-1
  }
  ninter<-max(down, ng-up)
  nMin<-ceiling(trim*ng)

  possibleThresh <- abs(allTh-threshRef)
  wh.thresh <- max(which(possibleThresh==min(possibleThresh)))
  
#search for a second threshold smaller than the first one
  if(wh.thresh>down+nMin){
    upInter<-wh.thresh-nMin
    if(allTh[upInter]==allTh[upInter-1])
      upInter<-which(allTh==allTh[upInter])[1]-1
    gammaMinus<-unique(allTh[seq(from=down+1, to=upInter)])
    #if only one unique value in middle regime
     if(allThUniq[which(allThUniq==allTh[upInter])+1]==allTh[wh.thresh]){
       gammaMinus <- head(gammaMinus, -1)#cut last one
      }
    if(length(gammaMinus)>0)
      storeMinus <- mapply(fun,gam1=gammaMinus,gam2=threshRef, thDelay=delayRef, MoreArgs=More)
    else
      storeMinus <- NA
  }
  else
    storeMinus <- NA
  
	#search for a second threshold higher than the first
  if(wh.thresh<up-nMin){
    downInter <- wh.thresh+nMin
    if(FALSE){#allTh[downInter]==allTh[wh.thresh]){
      samesTh <-which(allTh==allTh[downInter])
      downInter <-samesTh[length(samesTh)]+nMin
    }    
    if(allTh[downInter]==allTh[downInter+1]){
      samesInter <-which(allTh==allTh[downInter])
      downInter <-samesInter[length(samesInter)]+1
    }
    gammaPlus<-unique(allTh[seq(from=downInter, to=up)])
          #if only one unique value in middle regime
    if(allThUniq[which(allThUniq==allTh[downInter])-1]==allTh[wh.thresh]){
      gammaPlus <- gammaPlus[-1]#cut first one
    }
    if(length(gammaPlus)>1)
      storePlus <- mapply(fun,gam1=threshRef,gam2=gammaPlus, thDelay=delayRef,MoreArgs=More)
    else
      storePlus <- NA
  }
  else
    storePlus <- NA

	#results
  store2 <- c(storeMinus, storePlus)
  positionSecond <- which(store2==min(store2, na.rm=TRUE), arr.ind=TRUE)
  if(positionSecond<=length(storeMinus))
    newThresh<-gammaMinus[positionSecond]
  else
    newThresh<-gammaPlus[positionSecond-length(storeMinus)]
  SSR<-min(store2, na.rm=TRUE)
  if(trace)
    cat("Second best: ",newThresh, " (conditionnal on th= ",threshRef, " and Delay= ", delayRef," ) \t SSR/AIC: ", SSR, "\n", sep="")
  list(threshRef=threshRef, newThresh=newThresh, SSR=SSR)
}


###Function for nthresh=2, Alternative step: grid around the points from first step
grid<-function(gammasUp, gammasDown, fun, trace=TRUE, method=c("for", "apply", "mapply"),...){
	store <- matrix(NA,ncol=length(gammasUp), nrow=length(gammasDown))
	method<-match.arg(method)
	if(method=="for"){
		#Grid search
		for(i in seq_len(length(gammasDown))){
			gam1 <- gammasDown[i]
			for(j in 1: length(gammasUp)){
				gam2 <- gammasUp[j]
				store[i,j] <- fun(gam1=gam1, gam2=gam2,...)
			}
		}

		#Finding the best result
		positionIter <- which(store==min(store, na.rm=TRUE), arr.ind=TRUE)
		rIter <- positionIter[1]
		cIter <- positionIter[2]
	
		gamma1Iter <- gammasDown[rIter]
		gamma2Iter <- gammasUp[cIter]

		bestThresh <- c(gamma1Iter, gamma2Iter)
	}
	else if(method=="apply"){
		grid<-expand.grid(gammasDown,gammasUp)
# 		fun(gam1=gam1, gam2=gam2, d=bestDelay)
		temp<-function(a) fun(gam1=c(a)[1],gam2=c(a)[2],...)
		store<-apply(grid,1,temp)
		bests<-which(store==min(store, na.rm=TRUE))
		if(length(bests)>1) {
			warning("There were ",length(bests), " thresholds values which minimize the SSR in the first search, the first one was taken") 	
			bests<-bests[1]}
# 		beta_grid<-grid[bests,1]
# 		bestGamma1<-grid[bests,2]
		bestThresh <- c(grid[bests,1], grid[bests,2])
	}
	else if(method=="mapply"){
		grid<-expand.grid(gammasDown,gammasUp)	
		store<-mapply(fun, gam1=grid[,1],gam2=grid[,2], MoreArgs=list(...))
		bests<-which(store==min(store, na.rm=TRUE))
		if(length(bests)>1) {
			warning("There were ",length(bests), " thresholds values which minimize the SSR in the first search, the first one was taken") 	
			bests<-bests[1]}
# 		beta_grid<-grid[bests,1]
# 		bestGamma1<-grid[bests,2]
		bestThresh <- c(grid[bests,1], grid[bests,2])
	}
	if(trace)
		cat("\nSecond step best thresholds", bestThresh, "\t\t SSR:", min(store, na.rm=TRUE), "\n")
	return(bestThresh)
}#end of grid function
#  VAR<-TVAR(dat[1:300,], lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="apply", max.iter=5)

#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="apply"))
#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="mapply"))
#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="for"))

if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
#data(zeroyld)
dat<-zeroyld
environment(TVAR)<-environment(star)
environment(summary.TVAR)<-environment(star)

VAR<-TVAR(as.matrix(dat[1:100,]),lag=2, nthresh=1,thDelay=1,trim=0.1, plot=FALSE, commonInter=TRUE, include="const", gamma=c(3.016, 3.307))

VAR<-TVAR(dat[1:100,],lag=2, nthresh=2,thDelay=1,trim=0.1, plot=FALSE, commonInter=TRUE, include="const", gamma=c(3.016, 3.307))
VAR<-TVAR(dat[1:100,], lag=2, nthresh=2,thDelay=1,trim=0.1, plot=FALSE, commonInter=TRUE, include="const")
#lag2, 2 thresh, trim00.05: 561.46

class(VAR)
VAR
print(VAR)
logLik(VAR)
AIC(VAR)
BIC(VAR)
coef(VAR)
deviance(VAR)
summary(VAR)
toLatex(VAR)
toLatex(summary(VAR))
VAR[["model.specific"]][["oneMatrix"]]
###TODO
#pre specified gamma: not working!
#Bname: give different output... either matrix or list... to improve
}
