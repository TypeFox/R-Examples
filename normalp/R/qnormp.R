qnormp<-function(pr,mu=0,sigmap=1,p=2,lower.tail=TRUE,log.pr=FALSE){
if(!is.numeric(pr)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
if (log.pr==TRUE) pr<-log(pr)
if (lower.tail==FALSE) pr<-1-pr
zp<-ifelse(pr<0.5,0.5-pr,pr-0.5)
zp<-2*zp
qg<-qgamma(zp,shape=1/p,scale=p)
z<-qg^(1/p)
z<-ifelse(pr<0.5,-z,z)
q<-mu+z*sigmap
q
}

