#'Test of unit root against SETAR alternative
#'
#'Test of unit root against a stationnary three regime SETAR alternative
#'
#'TODO
#'
#'@param x time series
#'@param m Number of lags under the alternative
#'@param series time series name (optional)
#'@param testStat Type of test statistic to use
#'@param trim trimming parameter indicating the minimal percentage of
#'observations in each regime
#'@param grid Whether a minimal number of percentage or observations should be
#'imposed. See details
#'@return A object of class "BBC2004Test" containing:
#'
#'-The value of the sup Test
#'
#'-The version of test used (either Wald, LM or LR).
#'@author Matthieu Stigler
#'@seealso \code{\link{setarTest}} for a test with stationarity as a null.
#'@export
#'@examples
#'
#'BBCTest(lynx, m=3, test="Wald", grid="minPerc")
#'
BBCTest<-function(x, m, series, testStat=c("LR", "Wald", "LM"), trim=0.1, grid=c("minPerc", "minObs")){

warning("This function should be considered as experimental. No consistency checks with original paper could be made\n")
test<-match.arg(testStat)
grid<-match.arg(grid)
#demean the data
if(abs(mean(x))>0.01)
  x<-x-mean(x)

  if(missing(series))
    series <- deparse(substitute(x))
    
  str <- nlar.struct(x=x, m=m, d=1, steps=1, series=series)
  z<-matrix(getdX1(str), ncol=1)

### selectSETAR 3: set-up of the grid


nobs<-length(x)
down<-trim*nobs
up<-(1-trim)*nobs

thValues<-matrix(sort(abs(z))[seq(down, up)], ncol=1)


#check that there is at least enough observations in each regime for the regression
  if(trim*length(z)<m+2){
    warning("It may have not enough observations with this trim parameter, has been changed from ", trim, " to ", round((m+2)/length(z),3))
    trim<-(m+2)/length(z)
  }
  critMin<-switch(grid, "minPerc"=trim, "minObs"= m+2)
  fun<-switch(grid, "minPerc"=mean, "minObs"= sum)
  dep<-switch(grid, "minPerc"=1, "minObs"= length(z))
  A<-apply(thValues, 1, PercInRegime, Xtrans=z,fun=fun, crit=critMin,dep=dep)
  B<-cbind(thValues,A)
  
  th<-B[which(B[,2]==1),1]
  if(length(th)==0)
    stop("Not enough observations or trimming parameter too high")
###test
  LR<-function(gam){
    gam1<- -abs(gam)
    gam2<- abs(gam)
    yy <- getdYY(str)
    
    ###unrestricted model
    xxUnr <- cbind(getdX1(str),getdXX(str))
    SSRUnrest<-SSR_2threshNoCommon(gam1,gam2,thDelay=0, yy=yy,xx=xxUnr,trans=z, ML=1:(m+1), MH=1:(m+1), MM=1:(m+1),const=1,trim=trim)
  
    ###restricted model
    xxR <- getdXX(str)
    SSRRest<-SSR_2threshNoCommon(gam1,gam2,thDelay=0, yy=yy,xx=xxR,trans=z, ML=1:m, MH=1:m, MM=1:m,const=1,trim=trim)
  
    ###stat  
    stat<-nrow(xxR)*log(SSRRest/SSRUnrest)
    n<-nrow(xxR)
  
  
    return(stat)
  }

  Wald<-function(gam){
    gam1<- -abs(gam)
    gam2<- abs(gam)
    yy <- getdYY(str)
    
    #unrestricted
    xx <- cbind(getdX1(str),getdXX(str))
    XX<-buildXth2NoCommon(gam1,gam2,thDelay=0,xx=xx,trans=z, ML=1:(m+1), MH=1:(m+1), MM=1:(m+1),const=1,trim=trim)
    
    #regress, extract coef and var of resid
    lm<-lm.fit(XX,yy)
    co<-lm$coefficients
    resVar<-crossprod(lm$residuals)/nrow(xx)
    
    #matrix of restriction
    R<-matrix(0, nrow=3, ncol=3*m+6)
    R[1,2]<-1
    R[2,2+m+2]<-1
    R[3,2+m+2+m+2]<-1
    XXInv<-solve(crossprod(XX))
    Rco<-R%*%co
    
    #test stat
    stat<-t(Rco)%*%solve(R%*%XXInv%*%t(R))%*%Rco/resVar
    return(stat)
  
  }


LM<-function(gam){
  gam1<- -abs(gam)
  gam2<- abs(gam)
  yy <- getdYY(str)
  
  #unrestricted
  xx <- cbind(getdX1(str),getdXX(str))
  XX<-buildXth2NoCommon(gam1,gam2,thDelay=0,xx=xx,trans=z, ML=1:(m+1), MH=1:(m+1), MM=1:(m+1),const=1,trim=trim)
  XXInv<-solve(crossprod(XX))
  
  #restricted
  xxRes <- getdXX(str)
  xxResTh<-buildXth2NoCommon(gam1,gam2,thDelay=0, xx=xxRes,trans=z, ML=1:m, MH=1:m, MM=1:m,const=1,trim=trim)
  lm<-lm.fit(xxResTh,yy)
  residRes<-lm$residuals
  resVar<-crossprod(residRes)/nrow(xx)
  
  stat<-t(residRes)%*%XX%*%XXInv%*%t(XX)%*%residRes/resVar
  return(stat)
}
 

allTest<-apply(matrix(th, ncol=1), 1, test)
maxStat<-max(allTest)

ret<-list(statistic=maxStat, method=test)
class(ret)<-"BBC2004Test"
 return(ret)
}

print.BBC2004Test<-function(x, ...){
  cat("Test of unit root against stationary setar\n\n")
  
  #names(x$maxStat)<-paste("max", x$test, sep="")
  cat("Test statistic:") 
  metName<-paste("max", x$method, sep="")
  print(matrix(x$statistic, nrow=1, dimnames=list(metName,"")))
  cat("\nCritical values:\n")
  critW<-c(16.181, 18.400, 23.010)
  critLM<-c(15.587, 17.630, 21.756)
  critLR<-c(15.772, 17.898, 22.232)
  Crit<-switch(x$method, "LR"=critLR, "LM"=critLM, "Wald"=critW)
  print(matrix(Crit, nrow=1, dimnames=list("",c(0.9, 0.95,0.99))))

}



if(FALSE){
library(tsDyn)
environment(BBCTest)<-environment(star)


BBCTest(lynx, m=3, test="Wald", grid="minPerc")
RW<-cumsum(rnorm(100))
#BBCtest(RW, m=3, test="Wald")
##BBCtest(RW, m=3, test="LR")
#BBCtest(RW, m=3, test="LM")
}



if(FALSE){
  #problems:
library(tsDyn)
library(zoo)
a<-load(file="/home/mat/NIPFP/Applied/Arbitrage/adrdata.rda")
ratio<-icici$Ratio
ratioIc <- as.vector(ratio)
environment(BBCTest)<-environment(star)
  BBCTest(ratioIc, m=1)
  BBCTest(ratioIc, m=1, testStat="Wald")
  BBCTest(ratioIc, m=1, testStat="LM")
#negative value...


ar<-cumsum(rnorm(20))
  Wald <- BBCTest(ar, m=1, testStat="Wald")$statistic

}
PercInRegime<-function(th, Xtrans, fun, crit, dep){
  LowRegime<-fun(Xtrans<= -abs(th))
  HighRegime<-fun(Xtrans> abs(th))
  MiddleRegime<-dep-LowRegime-HighRegime
  RegProportions<-c(LowRegime, MiddleRegime, HighRegime)
  
  ret<-ifelse(any(RegProportions<crit), 0,1)
  return(ret)
}
