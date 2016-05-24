jackknifeKME <-
function(X, Y, delta, method="PDQ", estimator = 1)
{

km.est<-function(Y, delta, estimator)
  {
	 Ynew<-Y^estimator
	 sorted<-order(Y) 
  	 sdelta<-as.integer(delta[sorted])
  	 delta.nminus1<-sdelta[(length(sdelta)-1)]
   	 if(delta.nminus1==0)
	 W<-kmweight.corr(Y, delta)
      else
      W<-kmweight(Y, delta)
	 kme<-sum(W*Ynew)
	 return(kme)
  }


modkm.est<-function(X, Y, delta, method="PDQ", estimator)
  {
        n<-length(delta)
 	   sorted<-order(Y) 
  	   sdelta<-as.integer(delta[sorted])
  	   delta.nminus1<-sdelta[(n-1)]
   if(method=="PDQ")
     {phiYn.mod<-imputeYn(X,Y,delta,method="PDQ")$Yn
      sdelta[n]<-1
     }
   else if(method=="condMean")
     {phiYn.mod<-imputeYn(X,Y,delta,method="condMean")$Yn
      sdelta<-imputeYn(X,Y,delta, method="condMean")$newdata[,2]
     }    
   else if(method=="condMedian")
     {phiYn.mod<-imputeYn(X,Y,delta,method="condMedian")$Yn
      sdelta<-imputeYn(X,Y,delta, method="condMedian")$newdata[,2]
     }    
   else if(method=="RcondMean")
     {phiYn.mod<-imputeYn(X,Y,delta,method="RcondMean")$Yn
      sdelta<-imputeYn(X,Y,delta, method="RcondMean")$newdata[,2]
     }    
   else
     {phiYn.mod<-imputeYn(X,Y,delta,method="RcondMedian")$Yn
      sdelta<-imputeYn(X,Y,delta, method="RcondMedian")$newdata[,2]
     }    

	 Ynew<-Y^estimator
   	 if(delta.nminus1==0)
	 W<-kmweight.corr(Y, delta)
      else
      W<-kmweight(Y, delta)
	 kme<-sum(W*Ynew)
	 return(kme)
  }


phiYn<-function(Y)
  {
  sorted<-order(Y) 
  sY<-as.double(Y[sorted])
  phiY<-sY[length(sY)]   
return(phiY)
 }

Jbias.kme<-function(Y, delta, estimator)
  {
   srt<-order(Y)
   sy<-as.double(Y[srt])
   sdelta<-as.integer(delta[srt])
   n <- length(sdelta)
   if(n != length(sdelta) || n != length(Y))
   stop("dimensions of Y and delta don't match!")
   phiYn <- phiYn(Y)
   biasf <- numeric(n-2)
   biasf[1]<-((n-2)/(n-1))^sdelta[1]
   for(i in 2 : (n-2))
   {
   biasf[i]<- biasf[i-1]*(((n-i-1)/(n-i))^sdelta[i])
   }
   bias<--(n-1)/n*(phiYn^estimator)*sdelta[n]*(1-sdelta[n-1])*biasf[n-2]
   return(bias)
 }

Bcorr.Jkme<-function(Y, delta, estimator)
  {
    jkme<-km.est(Y, delta, estimator)-Jbias.kme(Y, delta, estimator)
    return(jkme) 
  }

modJbias.kme<-function(X, Y, delta, method="PDQ", estimator)
  {
   srt<-order(Y)
   sy<-as.double(Y[srt])
   sdelta<-as.integer(delta[srt])
   n <- length(sdelta)
   if(n != length(sdelta) || n != length(Y))
   stop("dimensions of Y and delta don't match!")
   if(method=="PDQ")
     {phiYn.mod<-imputeYn(X,Y,delta,method="PDQ")$Yn
      sdelta[n]<-1
     }
   else if(method=="condMean")
     {phiYn.mod<-imputeYn(X,Y,delta,method="condMean")$Yn
      sdelta<-imputeYn(X,Y,delta, method="condMean")$newdata[,2]
     }    
   else if(method=="condMedian")
     {phiYn.mod<-imputeYn(X,Y,delta,method="condMedian")$Yn
      sdelta<-imputeYn(X,Y,delta, method="condMedian")$newdata[,2]
     }    
   else if(method=="RcondMean")
     {phiYn.mod<-imputeYn(X,Y,delta,method="RcondMean")$Yn
      sdelta<-imputeYn(X,Y,delta, method="RcondMean")$newdata[,2]
     }    
   else
     {phiYn.mod<-imputeYn(X,Y,delta,method="RcondMedian")$Yn
      sdelta<-imputeYn(X,Y,delta, method="RcondMedian")$newdata[,2]
     }    
   biasf <- numeric(n-2)
   biasf[1]<-((n-2)/(n-1))^sdelta[1]
   for(i in 2 : (n-2))
   {
   biasf[i]<- biasf[i-1]*(((n-i-1)/(n-i))^sdelta[i])
   }
   bias<--(n-1)/n*(phiYn.mod^estimator)*sdelta[n]*(1-sdelta[n-1])*biasf[n-2]
   return(bias)
 }

Bcorr.modJkme<-function(X, Y, delta, method="PDQ", estimator)
  {
    modjkme<-modkm.est(X, Y, delta, method="PDQ", estimator)-modJbias.kme(X, Y, delta, method="PDQ", estimator)
    return(modjkme)
  }

jackBias=list(km.est = km.est(Y, delta, estimator), modkm.est= modkm.est(X, Y, delta, method = "PDQ", estimator), Jbias.kme = Jbias.kme(Y, delta, estimator), Bcorr.Jkme = Bcorr.Jkme(Y, delta, estimator), modJbias.kme = modJbias.kme(X, Y, delta, method="PDQ", estimator), Bcorr.modJkme = Bcorr.modJkme(X, Y, delta, method="PDQ", estimator))
class(jackBias)<-("jacknifeKME")
return(jackBias)
}
