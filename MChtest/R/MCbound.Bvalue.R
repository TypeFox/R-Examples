"MCbound.Bvalue" <-
function(nmax,alpha,e0,e1,conf.level=.99){

    theta<-findB.design(nmax,alpha,e0,e1)$theta
    theta<-sorttheta(theta)$theta
    out<-find.cb(theta)
    r<-out$w
    n<-out$m
    cb<-out$cb
    S<-length(r)
    cl.alpha<-1-conf.level
    check<-  sum( cb*mult(n,r,p=.03828483) )
    ci.lower<-ci.upper<-p.value<-rep(0,S)
    cl.alpha<- (1-conf.level)/2
    phat<-r/n
    #phat<-(r+1)/(n+1)
    eps<-.Machine$double.eps*2
    for (i in 1:S){
        ### first calculate p-value

        index<- phat <= phat[i]
        p.value[i]<-sum( cb[index] )
  
        func1<-function(p,I=i){
            index<- phat >= phat[I]
            np<-length(p)
            out<-rep(NA,np)
            for (j in 1:np){
                 out[j]<-sum( cb[index]*mult(n[index],r[index],p[j]) )-cl.alpha
            }
            out
        }
        if (r[i]==0){ ci.lower[i]<-0 }
        else { ci.lower[i]<-uniroot(func1,c(eps,1-eps))$root }
        func2<-function(p,I=i){
            index<- phat <= phat[I]
            np<-length(p)
            out<-rep(NA,np)
            for (j in 1:np){
                 out[j]<-sum( cb[index]*mult(n[index],r[index],p[j]) )-cl.alpha
            }
            out
        }
        if (r[i]==n[i]) ci.upper[1]<-1
        else { ci.upper[i]<-uniroot(func2,c(eps,1-eps))$root }
    }
    parms<-c(nmax,alpha,e0,e1)
    names(parms)<-c("nmax","alpha","e0","e1")
    out<-list(S=r,N=n,p.value=p.value,ci.lower=ci.lower,
       ci.upper=ci.upper,Kstar=cb,conf.level=conf.level,
       type="Bvalue",
       parms=parms,check=check)
    class(out)<-"MCbound"
    out
}

