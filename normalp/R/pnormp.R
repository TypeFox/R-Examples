pnormp<-function(q,mu=0,sigmap=1,p=2,lower.tail=TRUE,log.pr=FALSE){
if(!is.numeric(q)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
z<-(q-mu)/sigmap
zz<-abs(z)^p
zp<-pgamma(zz,shape=1/p,scale=p)
zp<-zp/2
zp<-ifelse(z<0,0.5-zp,0.5+zp)
if (lower.tail==FALSE) zp<-1-zp
if (log.pr==TRUE) zp<-log(zp)
zp
}

