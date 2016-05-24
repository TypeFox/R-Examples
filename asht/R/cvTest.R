cvTest <-
function(x,nullCV=1, 
    alternative=c("two.sided","less","greater"),conf.level=0.95,
    distn=c("normal","lognormal"),CVmax=10^6){

    distn<-match.arg(distn)
    alternative<-match.arg(alternative)
    mu<-mean(x)
    s<- sd(x)
    n<-length(x)
    dname<-deparse(substitute(x))

    if (distn=="lognormal"){
        if (any(x<=0)) stop("some x[i]<=0 so lognormal distributional assumption violated")
        trafo<-function(v){ sqrt( exp(v)-1)  } 
        invt<-function(tau){ log(tau^2+1) }
        out<-var1Test(log(x),nullVar=invt(nullCV),alternative=alternative,conf.level=conf.level)
        out$estimate<- trafo(out$estimate)
        names(out$estimate)<-"Estimate Coef of Var from Log-normal Model"
        out$conf.int<- trafo(out$conf.int)
        out$method<-"Log-normal test of Coefficient of Variation"
        out$null.value<- trafo(out$null.value)
    } else if (distn=="normal"){
        # Under normality assumption
        # Tn ~ non-central t with df=n-1 and ncp= sqrt(n)/tau
        #  where tau=population coefficient of variance 
        tn<- sqrt(n)*mu/s
        NCP<- sqrt(n)/nullCV
        # since tn= sqrt(n)/tauhat
        #  reject null: tau<=tauB for alternative: tau> tauB
        # when tn is small 
        pAG<- pt(tn,n-1,ncp=NCP)
        # reject  null: tau>=tauB for alternative: tau< tauB
        # when tn is small
        pAL<- pt(tn,n-1,ncp=NCP,lower.tail=FALSE)
        alpha<-1-conf.level
        if (alternative=="two.sided"){ 
            alpha<-alpha/2 
            p.value<-2*min(pAL,pAG)
        } else if (alternative=="less"){
            p.value<-pAL
        } else if (alternative=="greater"){
            p.value<-pAG
        }

        rootUpper<-function(CV){
            alpha - pt(tn,n-1,ncp=sqrt(n)/CV,lower.tail=FALSE)
        }
        rootLower<-function(CV){
            alpha - pt(tn,n-1,ncp=sqrt(n)/CV)
        }
        CV<- s/mu
        lower<-0
        upper<-Inf
        if (!(CV<0 | CV==Inf)){
            # find root unless sample CV<0 or =Inf
            # otherwise keep at Inf
            upper<-uniroot(rootUpper,c(CV,10^6))$root
        }
        if (CV!=0){
            lower<-uniroot(rootLower,c(0,CV))$root
        }
        ci<-c(lower,upper)
        attr(ci,"conf.level")<-conf.level
        names(CV)<-"coefficient of variation"

        out<-list(p.value=p.value,
          conf.int=ci,
          estimate=CV,
          null.value=nullCV,
          alternative=alternative,
          method="Normal Test of Coefficient of Variation")
    }
    out$statistic<-mu
    names(out$statistic)<-"mean"
    out$parameter<-s
    names(out$parameter)<-"standard deviation"
    out$data.name<-dname
    class(out)<-"htest"
    return(out)
}
