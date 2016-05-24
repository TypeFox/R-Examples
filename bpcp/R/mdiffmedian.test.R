mdiffmedian.test <-
function(x1,x2,nulldiff=0,
    alternative=c("two.sided","less","greater"),
    conf.level=0.95){
    beta0<-nulldiff
    alternative<-match.arg(alternative)

    n1<-length(x1)
    tab1<-table(x1)
    k1<-length(tab1)
    ux1<-sort(unique(x1))
    
    n2<-length(x2)
    tab2<-table(x2)
    k2<-length(tab2)
    ux2<-sort(unique(x2))

    d1<-cumsum(tab1)
    F1<-pbinom(c(0,d1),n1,.5)
    f1<-F1 - c(0,F1[-(k1+1)]) 
    L1<-c(-Inf,ux1)
    U1<-c(ux1,Inf)

    d2<-cumsum(tab2)
    F2<-pbinom(c(0,d2),n2,.5)
    f2<-F2 - c(0,F2[-(k2+1)]) 
    L2<-c(-Inf,ux2)
    U2<-c(ux2,Inf)

    probs<- matrix(f1,k1+1,1) %*% matrix(f2,1,k2+1)

    LL1<-matrix(rep(L1,k2+1),k1+1,k2+1)
    UU1<-matrix(rep(U1,k2+1),k1+1,k2+1)
    a<- (1-conf.level)/2
    lower1<- min(L1[F1>a])
    upper1<-min(U1[F1>1-a])

    LL2<-matrix(rep(L2,k1+1),k1+1,k2+1,byrow=TRUE)
    UU2<-matrix(rep(U2,k1+1),k1+1,k2+1,byrow=TRUE)

    lower2<- min(L2[F2>a])
    upper2<-min(U2[F2>1-a])

    LD<- as.vector(LL2 - UU1)
    UD<- as.vector(UU2 - LL1)
    pD<- as.vector(probs)
    oL<- order(LD)
    LD<-LD[oL]
    pL<-cumsum(pD[oL])
    lower<-min(LD[pL>a])

    oU<-order(UD)
    UD<-UD[oU]
    pU<-cumsum(pD[oU])
    upper<- min(UD[pU>1-a])
    
    pDL<-pD[oL]
    pvalL<- sum(pDL[LD<=beta0])
    pDU<-pD[oU]
    pvalU<- sum(pDU[UD>=beta0])
    pval2<- min(1,2*pvalL,2*pvalU)

    if (alternative=="less"){
        pvalue<-pvalU
    } else if (alternative=="greater"){
        pvalue<-pvalL
    } else if (alternative=="two.sided"){
        pvalue<-pval2
    }
    CI<-c(lower,upper)
    attr(CI,"conf.level")<-conf.level

    #out<-c(med1=median(x1),lower1=lower1,upper1=upper1,
    #        med2=median(x2),lower2=lower2,upper2=upper2,
    #        diff=median(x2)-median(x1),lower=lower,upper=upper,
    #        conf.level=conf.level,pvalL=pvalL,pvalU=pvalU,
    #   pval2=pval2)
    med1<-median(x1)
    names(med1)<-"median(x1)"
    med2<-median(x2)
    names(med2)<-"median(x2)"
    names(beta0)<-"median(x2)-median(x1)"
    dname<-paste("x1=",deparse(substitute(x1)),"and","x2=",deparse(substitute(x2)))
    est<-median(x2)-median(x1)
    names(est)<-"median(x2)-median(x1)"
    out<-list(
        statistic=med1,
        parameter=med2,
        p.value=pvalue,
        conf.int=CI,
        estimate=est,
        null.value=beta0,
        alternative=alternative,
        method="Melded Difference in Medians Test (using Sign Test CIs)",
        data.name=dname)
    class(out)<-"htest"
    #round(out,4)
    out
}
