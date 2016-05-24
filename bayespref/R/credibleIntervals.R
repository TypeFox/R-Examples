credibleIntervals <-
function(prefres=NULL,burn=0,interval=0.95){
  dif<-(1 - interval) / 2
  lb<-0 + dif
  ub<-1 - dif
  mcmc<-dim(prefres$PopPref)[2]
  nCat<-dim(prefres$PopPref)[1]
  nInd<-dim(prefres$IndPref)[1]
  ## ind pref
  if (is.null(prefres$IndPref)==FALSE){
    ip<-prefres$IndPref[,,(burn+1):mcmc]
    ## ind x catagories x lb med ub
    ipQ<-array(dim=c(nInd,nCat,3))
    for (i in 1:nInd){
      for (j in 1:nCat){
        ipQ[i,j,]<-quantile(ip[i,j,],probs=c(lb,0.5,ub))
      }
    }
  }
  else{
    ipQ<-NULL
  }
  ## pop pref
  pp<-prefres$PopPref[,(burn+1):mcmc]
  ppQ<-array(dim=c(nCat,3))
  for (j in 1:nCat){
    ppQ[j,]<-quantile(pp[j,],probs=c(lb,0.5,ub))
  }
  pm<-prefres$PopVar[(burn+1):mcmc]
  pmQ<-quantile(pm,probs=c(lb,0.5,ub)) 
  resout<-list(ipQ,ppQ,pmQ)
  names(resout)<-c("IndPref","PopPref","PopVar")
  return(resout)
}

