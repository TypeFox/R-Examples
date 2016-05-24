"MCtest.default" <-
function(x,statistic,resample,bound,extreme="geq",seed=1234325){
    set.seed(seed)
    t0<-statistic(x)
    minn<-min(bound$N)
    maxn<-max(bound$N)
    s<-0
    n<-0
    index<- 1:length(x)
    ti<-rep(NA,maxn)
    for (i in 1:minn){
        n<-n+1
        ti[n]<-statistic(resample(x))
    }
    if (extreme=="leq") s<-length((1:n)[ti[1:n]<=t0])
    else if (extreme=="geq") s<-length((1:n)[ti[1:n]>=t0])
    stop.index<-n==bound$N & s==bound$S
    if (any(stop.index)){
        if (length(bound$N[stop.index])==1){
            p<-bound$p.value[stop.index]
            ci.lower<-bound$ci.lower[stop.index]
            ci.upper<-bound$ci.upper[stop.index]
            pvalue.ci<-c(ci.lower,ci.upper)
            attr(pvalue.ci,"conf.level")<-bound$conf.level
        }
        else stop("stopping rule has duplicate values")
    }
    else{
        repeat{
            n<-n+1
            ti[n]<-statistic(resample(x))
            if (extreme=="leq") s<-ifelse(ti[n]<=t0,s+1,s+0)
            else if (extreme=="geq") s<-ifelse(ti[n]>=t0,s+1,s+0)
   
            stop.index<-n==bound$N & s==bound$S
            #cat("n=",n,"\n")
            #print(stop.index)
            if (any(stop.index) | n>maxn) break()
         }
         if (length(bound$N[stop.index])==1){ 
             p<-bound$p.value[stop.index]
             ci.lower<-bound$ci.lower[stop.index]
             ci.upper<-bound$ci.upper[stop.index]
             pvalue.ci<-c(ci.lower,ci.upper)
             attr(pvalue.ci,"conf.level")<-bound$conf.level
         }
         else if (length(bound$N[stop.index])>1) stop("stopping rule has duplicate values")
         else stop("reached n=max(bound$N) without stopping") 
    }
    out<-list(Ti=ti[!is.na(ti)],type=bound$type,parms=bound$parms,T0=t0,p.value=p,p.value.ci=pvalue.ci)
    class(out)<-"MCtest"
    out
}

