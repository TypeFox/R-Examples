#'Test of linearity
#'
#'Multivariate extension of the linearity against threshold test from Hansen
#'(1999) with bootstrap distribution
#'
#'This test is just the multivariate extension proposed by Lo and Zivot of the
#'linearity test of Hansen (1999). As in univariate case, estimation of the
#'first threshold parameter is made with CLS, for the second threshold a
#'conditional search with one iteration is made. Instead of a Ftest comparing
#'the SSR for the univariate case, a Likelihood Ratio (LR) test comparing the
#'covariance matrix of each model is computed.
#'
#'\deqn{ LR_{ij}=T( ln(\det \hat \Sigma_{i}) -ln(\det \hat \Sigma_{j}))} where
#'\eqn{ \hat \Sigma_{i}} is the estimated covariance matrix of the model with i
#'regimes (and so i-1 thresholds).
#'
#'Three test are avalaible. The both first can be seen as linearity test,
#'whereas the third can be seen as a specification test: once the 1vs2 or/and
#'1vs3 rejected the linearity and henceforth accepted the presence of a
#'threshold, is a model with one or two thresholds preferable?
#'
#'Test \bold{1vs}2: Linear VAR versus 1 threshold TVAR
#'
#'Test \bold{1vs}3: Linear VAR versus 2 threshold2 TVAR
#'
#'Test \bold{2vs3}: 1 threshold TAR versus 2 threshold2 TAR
#'
#'The both first are computed together and avalaible with test="1vs". The third
#'test is avalaible with test="2vs3".
#'
#'The homoskedastik bootstrap distribution is based on resampling the residuals
#'from H0 model, estimating the threshold parameter and then computing the
#'Ftest, so it involves many computations and is pretty slow.
#'
#'@aliases TVAR.LRtest TVAR.LRtest
#'@param data multivariate time series
#'@param lag Number of lags to include in each regime
#'@param trend whether a trend should be added
#'@param series name of the series
#'@param thDelay 'time delay' for the threshold variable (as multiple of
#'embedding time delay d) PLEASE NOTE that the notation is currently different
#'to univariate models in tsDyn. The left side variable is taken at time t, and
#'not t+1 as in univariate cases.
#'@param mTh combination of variables with same lag order for the transition
#'variable. Either a single value (indicating which variable to take) or a
#'combination
#'@param thVar external transition variable
#'@param nboot Number of bootstrap replications
#'@param plot Whether a plot showing the results of the grid search should be
#'printed
#'@param trim trimming parameter indicating the minimal percentage of
#'observations in each regime
#'@param test Type of usual and alternative hypothesis. See details
#'@param model Whether the threshold variable is taken in level (TAR) or
#'difference (MTAR)
#'@param hpc Possibility to run the bootstrap on parallel core. See details in
#'\code{\link{TVECM.HStest}}
#'@param trace should additional infos be printed? (logical)
#'@param check Possibility to check the function by no sampling: the test value
#'should be the same as in the original data
#'@return A list containing:
#'
#'-The values of each LR test
#'
#'-The bootstrap Pvalues and critical values for the test selected
#'@author Matthieu Stigler
#'@seealso \code{\link{setarTest}} for the univariate version.
#'\code{\link{OlsTVAR}} for estimation of the model.
#'@references Hansen (1999) Testing for linearity, Journal of Economic Surveys,
#'Volume 13, Number 5, December 1999 , pp. 551-576(26) avalaible at:
#'\url{http://www.ssc.wisc.edu/~bhansen/papers/cv.htm}
#'
#'Lo and Zivot (2001) "Threshold Cointegration and Nonlinear Adjustment to the
#'Law of One Price," Macroeconomic Dynamics, Cambridge University Press, vol.
#'5(4), pages 533-76, September.
#'@keywords ts
#'@export
#'@examples
#'
#'
#'data(zeroyld)
#'data<-zeroyld
#'
#'TVAR.LRtest(data, lag=2, mTh=1,thDelay=1:2, nboot=3, plot=FALSE, trim=0.1, test="1vs")
TVAR.LRtest <- function (data, lag=1, trend=TRUE, series, thDelay = 1:m, mTh=1, thVar, nboot=10, plot=FALSE, trim=0.1, test=c("1vs", "2vs3"), model=c("TAR", "MTAR"), hpc=c("none", "foreach"), trace=FALSE, check=FALSE) {

##Check args
  test<-match.arg(test)
  model<-match.arg(model)
  hpc<-match.arg(hpc)

  if (missing(series))  series <- deparse(substitute(data))

  if(is.null(colnames(data)))
	  colnames(data)<-paste("Var", seq_len(k), sep="")
  m<-lag			#keep consistent with nlar
  if(max(thDelay)>m)
	  stop("Max of thDelay should be smaller or equal to the number of lags")
  if(m<1)
	  stop("m should be at least 1")

  hasExternThVar <-ifelse(missing(thVar), FALSE, TRUE)

## Set-up variables
  y <- as.matrix(data) 
  ndig<-getndp(y)
  y<-round(y,ndig)
  Torigin <- nrow(y) 	#Size of original sample

  T <- nrow(y) 		#Size of start sample
  t <- T-lag 		#Size of end sample
  k <- ncol(y) 		#Number of variables
  p<-lag
  d<-1			#seems to be useless for nlVar...

  ndig<-getndp(y)
  Y <- y[(m+1):T,] #
  Z <- embed(y, m+1)[, -seq_len(k)]	#Lags matrix
  a<-0
  if(trend){
    Z <- cbind(1,Z)
    a<-1
  } else {
    warning("The test was currently implemented for a model with trend. Results could be altered without trend")
  }
  npar <- ncol(Z)		

  cat("Warning: the thDelay values do not correspond to the univariate implementation in tsdyn\n")

##################
###Linear model
#################
  B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar
  res<-Y-Z%*%t(B)
  Sigma<- matrix(1/T*crossprod(res),ncol=k)

  Y<-t(Y)

########################
### Threshold variable
########################

###External threshold variable
  if (hasExternThVar ) {		
	  if (length(thVar) != Torigin) {
		  z <- thVar[seq_len(Torigin)]
		  warning("The external threshold variable is not of same length as the original variable")
	  } else{
		  z <- thVar
	  }
	  z<-embed(z,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
  } else {

###Combination (or single value indicating position) of contemporaneous variables
    if (length(mTh) > k)
	    stop("length of 'mTh' should be equal to the number of variables, or just one")
    if(length(mTh)==1) {
	    if(mTh>k)
		    stop("mTh too big, should be smaller or equal to the number of variables")
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
	    z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]
  }
  z<-as.matrix(z)
  z<-round(z,ndig)


###############################
###Grid for transition variable
###############################

  if(length(thDelay)==1) 
	  b<-thDelay
  else 
	  b<-1
	  
  allgammas<-sort(z[,b])
  ng<-length(allgammas)
  # print(ng)
  nmin<-round(trim*ng)
  gammas<-unique(allgammas[ceiling(trim*ng+1):floor((1-trim)*ng-1)])


###################
###Search function
##################

  SSR_1thresh<- function(grid,Z,Y, trans){
	  d<-grid[1]
	  gam1<-grid[2]
	  ##Threshold dummies
	  d1<-ifelse(trans[,d]<gam1, 1,0)
	  ndown<-mean(d1)
	  if(min(ndown, 1-ndown)>=trim){
		  Z1 <- t(cbind(d1 * Z, (1-d1)*Z))		# dim k(p+1) x t
		  res<-crossprod(c(Y - tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))%*% Z1))}
	  else
		  res<-NA
	  return(res)
  } #end of the function

  Sigma_1thresh<- function(gam1, d,Z,Y, trans){
	  ##Threshold dummies
	  d1<-ifelse(trans[,d]<gam1, 1,0)
	  regimeDown <- d1 * Z
	  regimeUp<-(1-d1)*Z
	  ##SSR
	  Z1 <- t(cbind(regimeDown, regimeUp))		# dim k(p+1) x t
	  B1 <- tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))
	  matrix(1/T*tcrossprod(Y - B1 %*% Z1),ncol=k)
  } #end of the function

  SSR_2thresh <- function(gam1,gam2,d,Z,Y,trans){
	  ##Threshold dummies
	  dummydown <- ifelse(trans[,d]<gam1, 1, 0)
	  regimedown <- dummydown*Z
	  ndown <- mean(dummydown)
	  dummyup <- ifelse(trans[,d]>=gam2, 1, 0)
	  regimeup <- dummyup*Z
	  nup <- mean(dummyup)
	  ##SSR from TVAR(3)
  #  	print(c(ndown,1-nup-ndown,nup))
	  if(min(nup, ndown, 1-nup-ndown)>=trim){
		  Z2 <- t(cbind(regimedown, (1-dummydown-dummyup)*Z, regimeup))		# dim k(p+1) x t	
		  resid <- crossprod(c( Y - tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))%*%Z2))	#SSR
	  }
	  else
		  resid <- NA
	  return(resid)
  }

  Sigma_2thresh <- function(gam1,gam2,d,Z,Y,trans){
	  ##Threshold dummies
	  dummydown <- ifelse(trans[,d]<gam1, 1, 0)
	  regimedown <- dummydown*Z
	  dummyup <- ifelse(trans[,d]>=gam2, 1, 0)
	  regimeup <- dummyup*Z
	  ##SSR from TVAR(3)
	  Z2 <- t(cbind(regimedown, (1-dummydown-dummyup)*Z, regimeup))		# dim k(p+1) x t	
	  matrix(1/T*tcrossprod(Y - tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))%*%Z2),ncol=k)
  }

##################
###One threshold
##################

  IDS<-as.matrix(expand.grid(thDelay, gammas)) 
  result <- apply(IDS, 1, SSR_1thresh,Z=Z, Y=Y,trans=z)
  posBest<-which(result==min(result, na.rm=TRUE))
  bestDelay<-IDS[posBest,1]
  bestThresh<-IDS[posBest,2]

  if(trace) cat("Best unique threshold", bestThresh, "\t\t\t\t SSR", min(result, na.rm=TRUE), "\n")

  Sigma_mod1thresh<-Sigma_1thresh(gam1=bestThresh, d=bestDelay,Z=Z,Y=Y, trans=z)

##################
###Two thresholds
##################

###Function for conditional search
  condiStep<-function(allgammas, threshRef, delayRef, fun,MoreArgs=NULL, target=NULL){

    wh.thresh <- which.min(abs(allgammas-threshRef))
    Thr2<-which.min(abs(allgammas-target))
    
  #search for a second threshold smaller than the first
    if(wh.thresh>2*nmin){
      gammaMinus<-unique(allgammas[seq(from=max(nmin, Thr2-20), to=min(wh.thresh-nmin, Thr2+20))])
      storeMinus <- mapply(fun, gam1=gammaMinus,gam2=threshRef,MoreArgs=MoreArgs)	
    }
    else
      storeMinus <- NA
    
					  #search for a second threshold higher than the first
    if(wh.thresh<ng-2*nmin){
      gammaPlus<-unique(allgammas[seq(from=max(wh.thresh+nmin,Thr2-20), to=min(ng-nmin, Thr2+20))])
      storePlus <- mapply(fun,gam1=threshRef,gam2=gammaPlus,  MoreArgs=MoreArgs)
    }
    else
      storePlus <- NA
    
					  #results
    
    
    store2 <- c(storeMinus, storePlus)
    
    positionSecond <- which(store2==min(store2, na.rm=TRUE), arr.ind=TRUE)[1]
    if(positionSecond<=length(storeMinus))
      newThresh<-gammaMinus[positionSecond]
    else
      newThresh<-gammaPlus[positionSecond-length(storeMinus)]
    
  # cat("Second best: ",newThresh, " (conditionnal on ",threshRef, " ) \t SSR", min(store2, na.rm=TRUE), "\n")
    list(newThresh=newThresh, SSR=min(store2, na.rm=TRUE))
  }	#end function condistep


###Applying the function for conditional search to original data
  More<-list(d=bestDelay, Z=Z, Y=Y,trans=z)
  Thresh2<-condiStep(allgammas, bestThresh, fun=SSR_2thresh, MoreArgs=More)$newThresh
  Thresh3<-condiStep(allgammas, Thresh2, fun=SSR_2thresh, MoreArgs=More)
  smallThresh<-min(Thresh2, Thresh3$newThresh)
  bigThresh<-max(Thresh2, Thresh3$newThresh)

  if(trace){
    cat("Second best: ",Thresh2, " (conditionnal on ",bestThresh, ")\n")
    cat("Iterative best: ",Thresh3$newThresh, " (conditionnal on ",Thresh2, ")\n")
  }

  Sigma_mod2thresh<-Sigma_2thresh(gam1=smallThresh,gam2=bigThresh,d=bestDelay, Z=Z, Y=Y,trans=z)

###F test for original data
  LRtest12<-as.numeric(t*(log(det(Sigma))-log(det(Sigma_mod1thresh))))
  LRtest13<-as.numeric(t*(log(det(Sigma))-log(det(Sigma_mod2thresh))))
  LRtest23<-as.numeric(t*(log(det(Sigma_mod1thresh))-log(det(Sigma_mod2thresh))))
  LRs<-matrix(c(LRtest12, LRtest13, LRtest23),ncol=3, dimnames=list("Test", c("1vs2", "1vs3", "2vs3")))
  LRs <- if(test=="1vs") LRs[,-3]  else LRs[,3]

##############################
###Bootstrap for the F test
##############################

### Reconstruction of series from linear model for LRtest 1vs2 and 1vs3
  #initial data	
  res_lin<-res
  Yb<-matrix(0, nrow=nrow(y), ncol=k)		#Delta Y term
  Yb[1:m,]<-y[1:m,]			

  bootlinear<-function(res){ #res=res_lin
	  resi<-rbind(matrix(0,nrow=m, ncol=k),res[sample(seq_len(nrow(res)), replace=TRUE),])
	  if(check)
		  resi<-rbind(matrix(0,nrow=m, ncol=k),res)		#Uncomment this line to check the bootstrap
	  for(i in (m+1):(nrow(y))){
		  Yb[i,]<-rowSums(cbind(B[,1], B[,-1]%*%matrix(t(Yb[i-c(1:m),]), ncol=1),resi[i,]))
	  }
  return(Yb)
  }#end bootlinear
  # print(cbind(y, Yb))
### Reconstruction of series from 1 thresh model for LRtest 2vs3

  dummydown <- ifelse(z[,bestDelay]<=bestThresh, 1, 0)
  regimedown <- dummydown*Z
  Z2 <- t(cbind(regimedown, (1-dummydown)*Z))		# dim k(p+1) x t	
  B1thresh<-tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))
  res_thresh <- t(Y - B1thresh%*%Z2)	#SSR

  B1tDown<-B1thresh[,seq_len(k*m+1)]
  B1tUp<-B1thresh[,-seq_len(k*m+1)]

  #initial data	
  Yb2<-matrix(0, nrow=nrow(y), ncol=k)		#Delta Y term
  Yb2[1:m,]<-y[1:m,]			
  z2<-vector("numeric", length=nrow(y))
  z2[1:m] <- if(hasExternThVar) thVar[1:m] else y[1:m,]%*%combin	## fill first values of transition variable

  boot1thresh<-function(res){	#res=res_thresh
	  sampledNumbers<- sample(seq_len(nrow(res)), replace=TRUE)
	  resiT<-rbind(matrix(0,nrow=m, ncol=k),res[sampledNumbers,])
	  if(hasExternThVar) thVarBoot<-thVar[sampledNumbers] ## sample also the thVar!

	  if(check)
		  resiT<-rbind(matrix(0,nrow=m, ncol=k),res)

	  for(i in (m+1):(nrow(y))){
	    if(round(z2[i-bestDelay],ndig)<= bestThresh) {
	      Yb2[i,]<-rowSums(cbind(B1tDown[,1], B1tDown[,-1]%*%matrix(t(Yb2[i-c(1:m),]), 	ncol=1),resiT[i,]))
	    } else {
	      Yb2[i,]<-rowSums(cbind(B1tUp[,1], B1tUp[,-1]%*%matrix(t(Yb2[i-c(1:m),]), ncol=1),resiT[i,]))
	    }
	    z2[i]<-if(hasExternThVar) thVarBoot[i+(m+1)] else Yb2[i,]%*%combin
	  }
  # print(cbind(data,z2))
	  return(Yb2)
  }#end boot1thresh

  boot1threshBIS<-function(res){	
	  resiT<-rbind(matrix(0,nrow=m, ncol=k),res[sample(seq_len(nrow(res)), replace=TRUE),])
	  if(check)
		  resiT<-rbind(matrix(0,nrow=m, ncol=k),res)

	  for(i in (m+1):(nrow(y))){
		  if(round(z2[i-bestDelay]-z2[i-bestDelay-1],ndig)<= bestThresh) 
			  Yb2[i,]<-rowSums(cbind(B1tDown[,1], B1tDown[,-1]%*%matrix(t(Yb2[i-c(1:m),]), 	ncol=1),resiT[i,]))
		  else
			  Yb2[i,]<-rowSums(cbind(B1tUp[,1], B1tUp[,-1]%*%matrix(t(Yb2[i-c(1:m),]), ncol=1),resiT[i,]))
		  z2[i]<-(Yb2[i,])%*%combin
	  }
  return(Yb2)
  }#end boot1thresh

#####Bootstrap loop
  test<-switch(test, "1vs"="1vs", "2vs3"="2vs3")
  if(model=="TAR")
	  bootThresh<-boot1thresh
  else
	  bootThresh<-boot1threshBIS
	  
  bootModel<-switch(test, "1vs"=bootlinear, "2vs3"=bootThresh)
  resids<-switch(test, "1vs"=res_lin, "2vs3"=res_thresh)

  bootstraploop<-function(y){
    
    xboot<-round(bootModel(res=y),ndig)
    
    
  # Sigma of linear boot model
    string<-embed(xboot,m+1)
    Yboot <- string[,seq_len(k)] 	#
    Zb <- string[, -seq_len(k)]	#Lags matrix
    if(trend==TRUE)
      Zb <- cbind(1,Zb)
    Bboot<-t(Yboot)%*%Zb%*%solve(t(Zb)%*%Zb)		#B: OLS parameters, dim 2 x npar
    resboot<-Yboot-Zb%*%t(Bboot)
    Sigmab<- matrix(1/T*crossprod(resboot),ncol=k)
    
  #grid for threshold boot model
    
    if(hasExternThVar) {
      zcombin<-thVar
      zb<-embed(zcombin,p+1)[,seq_len(max(thDelay))+1]
    } else {
      zcombin<-  xboot%*%combin
      if(model=="MTAR"){
	if(max(thDelay)<p)
	  zb<-embed(diff(zcombin),p)[,seq_len(max(thDelay))+1]
	else if(max(thDelay)==p){
	  zb<-embed(diff(zcombin),p+1)[,seq_len(max(thDelay))+1]
	  zb<-rbind(0,as.matrix(zb))}
      }
      else
	zb <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]
    }
    zb<-as.matrix(zb)
    
  #  print(cbind(z,zb))

    allgammasb<-sort(zb[,b])
    gammasb<-unique(allgammasb[(ceiling(trim*ng)+1):floor((1-trim)*ng-1)])
    
  
###One threshold Search on bootstrap data

  IDSb<-as.matrix(expand.grid(thDelay, gammasb))
  resultb <- apply(IDSb, 1, SSR_1thresh,Z=Zb, Y=t(Yboot),trans=zb)
  postBestb<-which(resultb==min(resultb, na.rm=TRUE))[1]
  bestDelayb<-IDSb[postBestb,1]
  bestThreshb<-IDSb[postBestb,2]
  
  Sigma_mod1threshb<-Sigma_1thresh(gam1=bestThreshb, d=bestDelayb,Z=Zb,Y=t(Yboot), trans=zb)
# print(bestThreshb)
###Two threshold Search (conditional and 1 iteration) on bootstrap data
  Moreb<-list(d=bestDelayb, Z=Zb, Y=t(Yboot),trans=zb)
  Thresh2b<-condiStep(allgammasb, bestThreshb, fun=SSR_2thresh, MoreArgs=Moreb)$newThresh
  Thresh3b<-condiStep(allgammasb, Thresh2b, fun=SSR_2thresh, MoreArgs=Moreb)
  smallThreshb<-min(Thresh2b, Thresh3b$newThresh)
  bigThreshb<-max(Thresh2b, Thresh3b$newThresh)


  Sigma_mod2threshb<-Sigma_2thresh(gam1=smallThreshb,gam2=bigThreshb,d=bestDelayb, Z=Zb, Y=t(Yboot),trans=zb)
  
###Test statistic on bootstrap data
  LRtest12b<-as.numeric(t*(log(det(Sigmab))-log(det(Sigma_mod1threshb))))
  LRtest13b<-as.numeric(t*(log(det(Sigmab))-log(det(Sigma_mod2threshb))))
  LRtest23b<-as.numeric(t*(log(det(Sigma_mod1threshb))-log(det(Sigma_mod2threshb))))
  
  list(LRtest12b, LRtest13b, LRtest23b)
}#end of bootstraploop



  LRtestboot<-if(hpc=="none"){
    replicate(n=nboot,bootstraploop(y=resids))
  } else {
    foreach(i=1:nboot, .export="bootstraploop", .combine="rbind") %dopar% bootstraploop(y=resids)
  }

  LRtestboot12<-unlist(LRtestboot[1,])
  LRtestboot13<-unlist(LRtestboot[2,])
  LRtestboot23<-unlist(LRtestboot[3,])

  PvalBoot12<-mean(ifelse(LRtestboot12>LRtest12,1,0))
  CriticalValBoot12<-quantile(LRtestboot12, probs=c(0.9, 0.95, 0.975,0.99))
  PvalBoot13<-mean(ifelse(LRtestboot13>LRtest13,1,0))
  CriticalValBoot13<-quantile(LRtestboot13, probs=c(0.9, 0.95, 0.975,0.99))
  PvalBoot23<-mean(ifelse(LRtestboot23>LRtest23,1,0))
  CriticalValBoot23<-quantile(LRtestboot23, probs=c(0.9, 0.95, 0.975,0.99))

  if(test=="1vs"){
    CriticalValBoot<-rbind(CriticalValBoot12,CriticalValBoot13)
    PvalBoot<-c(PvalBoot12,PvalBoot13)
    names(PvalBoot)<-c("1vs2", "1vs3")
  } else {
    CriticalValBoot<-CriticalValBoot23
    PvalBoot<-PvalBoot23
  }


###Grahical output
#needs: LRtestboot12, LRtest12, m, k
  if(plot==TRUE&nboot>0){
	  if(test=="1vs"){
	    layout(c(1,2))
	    plot(density(LRtestboot12), xlab="LRtest12", xlim=c(0,max(LRtest12+1,max(LRtestboot12))),ylim=c(0,max(density(LRtestboot12)$y,dchisq(0:LRtest12, df=1+m))), main="Test linear VAR vs 1 threshold TVAR")
	    abline(v=LRtest12, lty=2, col=2)
	    dchis1 <- function(x) dchisq(x, df=1+k*m, ncp=0)
	    curve(dchis1, from=0, n=LRtest12+5, add=TRUE, col=3)
	    legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))

	    plot(density(LRtestboot13), xlab="LRtest13", xlim=c(0,max(LRtest13+1,max(LRtestboot13))),ylim=c(0,max(density(LRtestboot13)$y,dchisq(0:LRtest12, df=2*(1+m)))),main="Test linear VAR vs 2 thresholds TVAR")
	    abline(v=LRtest13, lty=2, col=2)
	    dchis2 <- function(x) dchisq(x, df=2*(1+k*m), ncp=0)
	    curve(dchis2,from=0, n=LRtest13+5, add=TRUE, col=3)
	    legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
	    layout(1)
	  }
	  else {
	    plot(density(LRtestboot23), xlab="LRtest23", xlim=c(0,max(LRtest23+1,LRtestboot23)), ylim=c(0,max(density(LRtestboot23)$y,dchisq(0:LRtest12, df=1+m))), main="Test 1 threshold TVAR vs 2 thresholds TVAR")
	    abline(v=LRtest23, lty=2, col=2)
	    dchis <- function(x) dchisq(x, df=1+k*m, ncp=0)
	    curve(dchis, from=0, n=LRtest23+5, add=TRUE, col=3)
	    legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
	  }
  }


#### RETURN RESULT
#nlar=extend(nlar(str, coef = res$coef, fit = res$fitted.values, res = res$residuals, k = res$k,
#list( model.specific = res),"setar")
  res<-list(bestDelay=bestDelay, LRtest.val=LRs, Pvalueboot=PvalBoot, CriticalValBoot=CriticalValBoot, type="test")
  class(res)<-"TVARtest"
  return(res)
}#End of thw whole function

print.TVARtest<-function(x,...){
  cat("Test of linear VAR against TVAR(1) and TVAR(2)\n\nLR test:\n")
  LR<-rbind(x$LRtest.val,x$Pvalueboot)
  rownames(LR)<-c("Test", "P-Val")
  print(LR)
}

summary.TVARtest<-function(object,...){
  cat("Test of linear VAR against TVAR(1) and TVAR(2)\n\nLR test:\n")
  LR<-rbind(object$LRtest.val,object$Pvalueboot)
  rownames(LR)<-c("Test", "P-Val")
  print(LR)
  cat("\n Bootstrap critical values for test 1 vs 2 regimes\n")
  print(object$CriticalValBoot[1,])
  cat("\n Bootstrap critical values for test 1 vs 3 regimes\n")
  print(object$CriticalValBoot[2,])
}


if(FALSE){ #usage example
environment(TVAR.LRtest)<-environment(star)
#data(zeroyld)
data<-zeroyld[1:150,]

test<-TVAR.LRtest(data, lag=3, mTh=c(1,1),thDelay=1:2, nboot=2, plot=FALSE, trim=0.1, test="1vs", model="TAR")
class(test)
print(test)
summary(test)
###Todo
environment(TVAR.LRtest)<-environment(star)
test<-TVAR.LRtest(data, lag=3, mTh=c(1,1),thDelay=1:2, nboot=2, thVar=data[,1], plot=FALSE, trim=0.1, test="1vs", model="TAR")


#does not work TVAR.LRtest(data, lag=3, mTh=c(1,1),thDelay=1:2, nboot=2, plot=TRUE, trim=0.1, test="1vs", check=TRUE, model="MTAR")
#
environment(TVAR.LRtest)<-environment(star)
TVAR.LRtest(data, lag=3, mTh=c(1,1),thDelay=1:2, nboot=2, plot=TRUE, trim=0.1, test="2vs3", check=TRUE, model="MTAR")
}


