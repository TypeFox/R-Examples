dnormp<-function(x,mu=0,sigmap=1,p=2,log=FALSE){
if(!is.numeric(x)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
cost<-2*p^(1/p)*gamma(1+1/p)*sigmap
expon1<-(abs(x-mu))^p
expon2<-p*sigmap^p
dsty<-(1/cost)*exp(-expon1/expon2)
if(log==TRUE) dsty<-log(dsty)
dsty
}

