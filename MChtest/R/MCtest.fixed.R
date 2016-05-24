"MCtest.fixed" <-
function(x,statistic,resample,Nmax=3619,extreme="geq",conf.level=.99,seed=1234325){
    N<-Nmax
    set.seed(seed)
    t0<-statistic(x)
    ti<-rep(NA,N)
    for (i in 1:N){
        ti[i]<-statistic(resample(x))
    }
    if (extreme=="leq") S<-length((1:N)[ti<=t0])
    else if (extreme=="geq") S<-length((1:N)[ti>=t0])

    p.value<-(S+1)/(N+1)
    p.alpha<-1-conf.level
    if (S==0){ ci.lower<-0 }
    else { ci.lower<- qbeta(p.alpha/2,S,N-S+1) }
    if (S==N){ ci.upper<-1 }
    else{ ci.upper<-qbeta(1-p.alpha/2,S+1,N-S) }

    pvalue.ci<-c(ci.lower,ci.upper)
    attr(pvalue.ci,"conf.level")<-conf.level

    out<-list(Ti=ti,type="fixed",parms=c(Nmax=N),T0=t0,p.value=p.value,p.value.ci=pvalue.ci)
    class(out)<-"MCtest"
    out
}

