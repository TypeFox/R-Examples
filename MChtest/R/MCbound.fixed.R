"MCbound.fixed" <-
function(Nmax,conf.level=.99){
    S<-0:Nmax
    N<-rep(Nmax,Nmax+1)
    p.value<- (S+1)/(N+1)
    alpha<-1-conf.level
    ci.lower<- c(0,qbeta(alpha/2,S[-1],N[-1]-S[-1]+1))
    ci.upper<- c(qbeta(1-alpha/2,S[-(Nmax+1)]+1,N[-(Nmax+1)]-S[-(Nmax+1)]),1)
    cb<-rep(1/(Nmax+1),Nmax+1)
    parms<-c(Nmax)
    names(parms)<-c("Nmax")
    out<-list(S=S,N=N,p.value=p.value,ci.lower=ci.lower,
       ci.upper=ci.upper,Kstar=cb,conf.level=conf.level,
       type="fixed",parms=parms)
    class(out)<-"MCbound"

    out
}

