exactbinomPlot<-function(x,
    n,
    p=NULL,
    ndiv=1000,
    tsmethod="central",
    pRange=c(0,1),
    dolines=FALSE,
    dopoints=TRUE,
    doci=TRUE, 
    alternative = c("two.sided", "less", "greater"), 
    relErr=1 + 10^(-7), 
    conf.level=.95,
    alphaline=TRUE,
    newplot=TRUE,
    midp=FALSE,...){

    ## copied setup code from binom.test, Version 2.10.1
    xr <- round(x)
    if (any(is.na(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
        stop("'x' must be nonnegative and integer")
    x <- xr
    if (length(x) == 2) {
        n <- sum(x)
        x <- x[1L]
    }
    else if (length(x) == 1) {
        nr <- round(n)
        if ((length(n) > 1) || is.na(n) || (n < 1) || abs(n - 
            nr) > 1e-07 || (x > nr)) 
            stop("'n' must be a positive integer >= 'x'")
        n <- nr
    }
    else stop("incorrect length of 'x'")


    if (is.null(p)){
        minp<-(pRange[1])
        maxp<-(pRange[2])
        p<-(minp + (0:ndiv)*(maxp-minp)/ndiv )
    }
    alternative <- match.arg(alternative)
    if (midp){
        pval <- switch(alternative, 
            less = pbinom(x, n, p)-0.5*dbinom(x,n,p), 
            greater = pbinom(x-1,n,p,lower.tail=FALSE)-
                0.5*dbinom(x,n,p), 
            two.sided = exactbinomPvals(x,n,p,relErr=relErr,
                tsmethod=tsmethod, midp=TRUE)$pvals)
    } else {
        pval <- switch(alternative, 
            less = pbinom(x, n, p), 
            greater = pbinom(x-1,n,p,lower.tail=FALSE), 
            two.sided = exactbinomPvals(x,n,p,relErr=relErr,
                tsmethod=tsmethod, midp=FALSE)$pvals)
    }


    if (newplot){
        YLAB<-"p-value"
        if (midp) YLAB<-"mid p-value"
        plot(p,pval,xlab="Null Hypothesis p",ylab=YLAB,
            type="n",...)
    } 
    if (dopoints){
        points(p,pval,...)
    } 
    if (dolines){
        lines(p,pval,...)
    }
    if (doci){
        ci<-binom.exact(x,n,tsmethod=tsmethod,alternative=alternative,
            control=binomControl(relErr=relErr),
            conf.level=conf.level,midp=midp)$conf.int
        alpha<-1-conf.level
        if (alphaline) lines(c(-1,2),c(alpha,alpha),lty=2)
        lines(c(ci[1],ci[1]),c(0,alpha),...)
        lines(c(ci[2],ci[2]),c(0,alpha),...)
    }
}

#exactbinomPlot(5,17,dolines=TRUE,dopoints=FALSE,ylim=c(0,.15))
#exactbinomPlot(5,17,tsmethod="blaker",col="blue",pch=".",doci=FALSE,newplot=FALSE)
#exactbinomPlot(5,17,tsmethod="blaker",col="blue",lty=2,dopoints=FALSE,newplot=FALSE)

