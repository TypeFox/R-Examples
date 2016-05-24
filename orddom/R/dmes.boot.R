### Calculate the SE and construct a CI for Cohen's d and the statistics used in dmes using bootstrap methods
### Based on Ruscio & Mullen's (2012) Bootstrap.SE.CI.A function (see Manual)
dmes.boot <- function(x,y,theta.es="dc",ci.meth="BCA",B=1999,alpha=.05,seed=1)
{  # initialize variables
   es  <- names(dmes(c(1:3),c(1:3)))
   ci  <- c("BSE","BP","BCA")
   if(is.element(theta.es,es) && is.element(ci.meth,ci)){
   theta<-match(theta.es,es)
   set.seed(seed)
   nx <- length(x)
   ny <- length(y)
   theta.obs <- unlist(dmes(x,y)[theta])
   d.obs <- metric_t(x,y,alpha)[5]
   CI.Lower <- CI.Upper <- theta.obs
   dCI.Lower<- dCI.Upper <- d.obs
   # perform bootstrap to generate B values of theta and t
   BS.Values <- rep(0,B)
    d.Values <- rep(0,B)
   for (i in 1:B) {
     xs<-sample(x,replace=TRUE)
	 ys<-sample(y,replace=TRUE)
     BS.Values[i]<-unlist(dmes(xs,ys)[theta])
	 d.Values[i] <-metric_t(xs,ys,alpha)[5]
	}
   BS.Values <- sort(BS.Values)
    d.Values <- sort(d.Values)
   # if all bootstrap samples yield same value for A, use it for both ends of CI
   if (min(d.Values) == max(d.Values))
      {dCI.Lower <- dCI.Upper <- d.Values[1]}
   if (min(BS.Values) == max(BS.Values))
      {CI.Lower <- CI.Upper <- BS.Values[1]} else {
   if (ci.meth=="BSE") {
   # construct trditional CI based on bootstrap SE
      CI.Lower <- theta.obs - (qnorm((1-(alpha/2))) * sd(BS.Values))
      CI.Upper <- theta.obs + (qnorm((1-(alpha/2))) * sd(BS.Values))
      dCI.Lower<- d.obs - (qnorm((1-(alpha/2))) * sd(d.Values))
      dCI.Upper<- d.obs + (qnorm((1-(alpha/2))) * sd(d.Values))
   }
   if (ci.meth=="BP") {	
   # bootstrap percentile CI
   # cf. Efron & Tibshirani (1993, Ch. 13) 
      CI.Lower <- BS.Values[round((alpha / 2) * B)]
      CI.Upper <- BS.Values[round((1 - alpha / 2) * B)]
      dCI.Lower <- d.Values[round((alpha / 2) * B)]
      dCI.Upper <- d.Values[round((1 - alpha / 2) * B)]
   }
   if (ci.meth=="BCA") {
   # bootstrap percentile bias corrected and accelerated
   # cf. Efron & Tibshirani (1993, Ch. 14) 
     # calculate bias-correction and acceleration parameters (z0 and a)
      z0 <- qnorm(mean(BS.Values < theta.obs)) # cf. Efron & Tibshirani (1993, p. 186, formula 14.14) 
      z0d<- qnorm(mean(d.Values < d.obs))
      jk <- rep(0, (nx + ny))
	  djk<- rep(0, (nx + ny))
      for (i in 1:nx)
         {jk[i] <- unlist(dmes(x[-i], y)[theta]) #jackknife
		  djk[i]<- metric_t(x[-1],y,alpha)[5]}
      for (i in 1:ny)
         {jk[nx+i] <- unlist(dmes(x, y[-i])[theta])
		  djk[nx+i]<- metric_t(x,y[-i],alpha)[5]}
      Diff <-mean(jk)-jk
	  Diffd<-mean(djk)-djk
      a<-sum(Diff^3)/(6*(sum(Diff^2))^1.5) # cf. Efron & Tibshirani (1993, p.186/15.15 and p. 328/22.29) 
	  ad<-sum(Diffd^3)/(6*(sum(Diffd^2))^1.5)
     # adjust location of endpoints, cf. Efron & Tibshirani (1993, p.185/14.10) 
      alpha1 <- pnorm(z0+((z0+qnorm(alpha/2))/(1-(a*(z0+qnorm(alpha/2))))))
      alpha2 <- pnorm(z0+((z0-qnorm(alpha/2))/(1-(a*(z0-qnorm(alpha/2))))))
      alpha1d<- pnorm(z0d+((z0d+qnorm(alpha/2))/(1-(ad*(z0d+qnorm(alpha/2))))))
      alpha2d<- pnorm(z0d+((z0d-qnorm(alpha/2))/(1-(ad*(z0d-qnorm(alpha/2))))))
      # if either endpoint undefined, replace it with value for percentile CI
      if (is.na(alpha1)) {alpha1<-(alpha/2)}
	  if (is.na(alpha1d)) {alpha1d<-(alpha/2)}
      if (is.na(alpha2)) {alpha2<-(1-(alpha/2))}
	  if (is.na(alpha2d)) {alpha2d<-(1-(alpha/2))}
      if (round(alpha1*B)<1) {CI.Lower <- BS.Values[1]} else {CI.Lower <- BS.Values[round(alpha1 * B)]}
	  if (round(alpha1d*B)<1) {dCI.Lower <- d.Values[1]} else {dCI.Lower <- d.Values[round(alpha1d*B)]}
      CI.Upper <- BS.Values[round(alpha2 * B)]	
	  dCI.Upper <- d.Values[round(alpha2d* B)]	
   }
   }
   # return A, SE of A, lower limit of CI, upper limit of CI
   return(list(theta=theta.obs,theta.SE=sd(BS.Values),bci.meth=ci.meth,theta.bci.lo=CI.Lower,theta.bci.up=CI.Upper,Coh.d=d.obs,Coh.d.bSE=sd(d.Values),Coh.d.bci.lo=dCI.Lower,Coh.d.bci.up=dCI.Upper))
}else{return(warning("theta.es or ci.meth parameter error"))}
}

