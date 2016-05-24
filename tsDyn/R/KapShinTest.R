
#' @export
KapShinTest <- function(x, m=1, series, include = c("none","const", "trend", "both"), c=3, delta=0.5, points=NULL,minObsMid=10, trick=c("for", "apply", "mapply"), trace=FALSE){
 
 warning("This function should be considered as experimental. No consistency checks with original paper could be made\n")
  include<-match.arg(include)
  if(missing(series))
    series <- deparse(substitute(x))

  
  #demean, demean detrend
  if(include!="none"){
    if(include=="trend")
      cat("Case of only detrending without demeaning is not considered in paper")
    const<-buildConstants(include=include, n=length(x))$const#in miscSETAR.R
    x <- residuals(lm(x~const))
  }


  
  #matrix of regressors    
  str <- nlar.struct(x=x, m=m, d=1, steps=1, series=series)
  xx <- cbind(getdX1(str),getdXX(str))
  yy <- getdYY(str)
  str$xx<-xx
  str$yy<-yy


  Y <- getdYY(str)
  XMinus1 <- getdX1(str)
  Xlags <- getdXX(str)
  
  

  R<-cbind(diag(1,2), matrix(0, nrow=2, ncol=m)) #restriction matrix

  bigY<-cbind(Y, XMinus1, Xlags)
  #Computation of wald for th1 and th2 given

#   wald <- function(gam1, gam2, Y=Y, Xminus1=Xminus1, Xlags=Xlags, const=const){
  wald <- function(gam1, gam2, Y, minObsMid){
    isL <- ifelse(Y[,2]<=gam1,1,0)
    isH <- ifelse(Y[,2]>gam2,1,0)
    if(nrow(Y)-sum(isL)-sum(isH)<minObsMid)
      return(NA)
    else{
      X <- cbind(isL*Y[,2], isH*Y[,2], Y[,-c(1:2)])
      reg <- lm.fit(X,Y[,1])
      
      beta <- matrix(reg$coefficients, ncol=1)    
      Rco<-R%*%beta
      SSR <- crossprod(reg$residuals)
      resVar <- SSR/length(beta)
      Wald <- (t(Rco)%*%solve(R%*%solve(crossprod(X))%*%t(R))%*%Rco)/resVar
      return(Wald)
    }
  }

# browser()
###grid
  dx1<-sort(getdX1(str))
  ndx1<-length(dx1)
  pbar<-0.5
  p1<-pbar-c/ndx1^delta
  p2<-pbar+c/ndx1^delta
  points<-if(is.null(points)) p2*ndx1-p1*ndx1 else points
  ths<-dx1[seq(from=p1*ndx1, to=p2*ndx1, length.out=points)]
  if(all(is.na(ths)))
    warning("Error in the grid specification")
  thMinus<-ths[1:round(0.5*length(ths))]
  thPlus<-ths[round(0.5*length(ths)):length(ths)]

###apply the grid
  GridSearch<-grid2(gammasUp=thPlus, gammasDown=thMinus, fun=wald, trace=trace, method="for",Y=bigY, minObsMid=minObsMid)
  walds<-GridSearch$store
  if(trace)
    cat(length(na.omit(walds)) ,"Combinations of thresholds were computed\n")
  SupW <- max(na.omit(walds))
  AvgW <- mean(na.omit(walds))
  ExpW <- mean(exp(na.omit(walds/2)))

  ret<-list(statistic=c(SupW, AvgW, ExpW), case=include, series=series)
  class(ret)<-"KapShin2006Test"
  return(ret)

}

print.KapShin2006Test<-function(x, ...){
  cat("Test of unit root against stationary setar\n\n")
  
  #names(x$maxStat)<-paste("max", x$test, sep="")
  cat("Test statistics:\n") 
  metName<-c("Sup","Ave", "ExpAve")
  print(matrix(x$statistic, nrow=1, dimnames=list("",metName)))
  cat("\nCritical values:\n")
   critVal<- matrix(c(6.01,7.29,10.35, 20.18, 38.28, 176.80,7.49,9.04, 12.16,42.30,91.83,  437.03,10.94, 12.64, 16.28, 237.46, 555.57, 3428.92), nrow=3, byrow=TRUE)
  Col<-switch(x$case, "none"=c(1,4), "const"=c(2,5), "both"=c(3,6))
  print(matrix(critVal[,Col], nrow=2, dimnames=list(c("Sup,Ave", "ExpAve"),c(0.9, 0.95,0.99))))
  cat("Data: ", x$series, switch(x$case, "none"="(raw)", "const"="(demeaned)", "both"="(demeaned and detrended)"), "\n")

}


if(FALSE){
library(tsDyn)
environment(KapShinTest) <- environment(star)
KS<-KapShinTest(lynx, m=1, trace=FALSE, include="none")
KS<-KapShinTest(lynx, m=1, trace=FALSE, points=10, include="none")
KS$statistic
KS
}

grid2<-function(gammasUp, gammasDown, fun, trace=TRUE, method=c("for", "apply", "mapply"),...){
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
	return(list(bestThresh=bestThresh, store=store))
}#end of grid function
#  VAR<-TVAR(dat[1:300,], lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="apply", max.iter=5)

#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="apply"))
#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="mapply"))
#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="for"))

  
