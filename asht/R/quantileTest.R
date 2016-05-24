quantileTest <- function(x,...){
    UseMethod("quantileTest")
}
quantileTest.default <-
function(x,q=0,prob=0.5,alternative=c("two.sided","less","greater"),
    conf.level=0.95,...){
    if (!(is.numeric(x) | is.integer(x))) stop("x must be numeric, integer, or an ordered factor")
    alternative <- match.arg(alternative)
    n<-length(x)
    ### p-values
    sLE<- length(x[x<=q])
    pAG<- pbinom(sLE,n,prob) 
    sGE<- length(x[x>=q])
    pAL<- pbinom(sGE,n,1-prob)
    pc<-min(1,2*pAG,2*pAL)
    if (alternative=="two.sided"){
        pval<-pc
    } else if (alternative=="less"){
        pval<-pAL
    } else if (alternative=="greater"){
        pval<-pAG
    }
    ## confidence intervals
    ux<-sort(unique(x))
    k<- length(ux)
    y<-rep(0,k+1)
    ystar<-rep(0,k+1)
    for (i in 1:k){
        y[i+1]<- length(x[x<=ux[i]])
        ystar[i]<-length(x[x>=ux[i]])
    }
    P<-pbinom(y,n,prob)   
    Pbar<-pbinom(ystar,n,1-prob)
    alpha<-1-conf.level
    if (alternative=="two.sided") alpha<-alpha/2
    # set extreme limits, in case we only calculate one-sided limits 
    lower<--Inf
    upper<- Inf
    if (alternative=="two.sided" | alternative=="greater"){
    # Lower Limits
        low<-c(-Inf,ux[1:k])
        lower<-min(low[alpha<P])
        #lower<-low[alpha<P & c(0,P[-(k+1)]) <= alpha]
        if (length(lower)!=1) stop("problem with lower limit code") 
    }
    if (alternative=="two.sided" | alternative=="less"){
    # Upper Limits
        hi<-c(ux[1:k],Inf)
        #lev<- 1-alpha
        #upper<-hi[lev<=P & c(0,P[-(k+1)]) < lev]
        upper<-max(hi[alpha<Pbar])
        if (length(upper)!=1) stop("problem with upper limit code") 
    }

    CI<-c(lower,upper)
    attr(CI,"conf.level")<- conf.level
    dname <- deparse(substitute(x))

    est<-quantile(x,probs=prob)
    names(q)<-names(est)<-ifelse(prob==0.5,"median","quantile")
    
    method<-ifelse(prob==0.5,
          "Exact Test/Confidence Interval for Median",
      "Exact Test/Confidence Interval for Quantile")

    statistic<- prob
    names(statistic)<-"quantile for prob"
    out<-list(statistic =statistic, parameter =c(pAG=pAG,pAL=pAL,pc=pc),
        p.value = pval, 
        conf.int = CI, estimate = est, null.value = q, 
        alternative = alternative, method = method, data.name = dname)
    class(out)<-"htest"
    out
}



quantileTest.ordered<-function(x,...){
    levx<-levels(x)
    xnum<-as.numeric(x)
    out<-quantileTest(xnum,...)
    getlevels<-function(x){
        nx<-length(x)
        xchar<-rep(NA,nx)
        floorx<-floor(x)
        ceilingx<-ceiling(x)
        xchar<- levx[floorx]
        inbtwn<- floorx!=ceilingx
        xchar[inbtwn]<- paste0(levx[floorx[inbtwn]],"/",levx[ceilingx[inbtwn]])
        attributes(xchar)<-attributes(x)
        xchar
     }

     out$conf.int<-getlevels(out$conf.int)
     out$estimate<-getlevels(out$estimate)
     out$null.value<-getlevels(out$null.value)
     out
}


medianTest <-
function(x,m=0,...){
    out<-quantileTest(x,q=m,prob=0.5,...)
    out
}



