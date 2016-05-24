exactpoissonPlot<-function(x,
    T=1,
    r=NULL,
    ndiv=1000,
    tsmethod="central",
    rRange=NULL,
    dolog=TRUE,
    dolines=FALSE,
    dopoints=TRUE,
    doci=TRUE, 
    alternative = c("two.sided", "less", "greater"), 
    relErr=1 + 10^(-7), 
    conf.level=.95,
    alphaline=TRUE,
    newplot=TRUE,
    midp=FALSE,...){

    ## copied setup code from poisson.test, Version 2.10.1
    ## copy much of this from poisson.test 
    ## for the two sample case primarily just change binom.test to binom.exact
    if ((l <- length(x)) != length(T)) 
        if (length(T) == 1) 
            T <- rep(T, l)
        else stop("'x' and 'T' have incompatible length")
    xr <- round(x)
    if (any(!is.finite(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
        stop("'x' must be finite, nonnegative, and integer")
    x <- xr
    if (any(is.na(T) | (T < 0))) 
        stop("'T' must be nonnegative")
    if ((k <- length(x)) < 1) 
        stop("not enough data")
    if (k > 2) 
        stop("The case k > 2 is unimplemented")
    alternative <- match.arg(alternative)


    ## 
    if (is.null(r) & is.null(rRange)){
        ## central method is the fastest method and most conservative
        ## we do not need to match the tsmethod
        ## we are just getting a rough range for the plot
        xtemp<-ifelse(x==0,1,x)
        ci<-poisson.exact(x=xtemp,T=T,r=1,
            tsmethod="central",alternative = "two.sided", 
            conf.level=.95,midp=midp)$conf.int
        rRange<-c(ci[1]*.8,ci[2]*1.2)
    }
    if (dolog & is.null(r)){
        minr<-log10(rRange[1])
        maxr<-log10(rRange[2])
        r<-10^(minr + (0:ndiv)*(maxr-minr)/ndiv )
    } else if (is.null(r)){
        minr<-(rRange[1])
        maxr<-(rRange[2])
        r<-(minr + (0:ndiv)*(maxr-minr)/ndiv )
    }


    if (k == 2) {
        p<- r* T[1]/(r * T[1] + T[2])
        n<-sum(x)
        midp.adjustment<- 0 
        if (midp) midp.adjustment<- 0.5*dbinom(x[1],n,p)
        if (midp & alternative=="two.sided" & tsmethod!="central") 
            stop("mid-p method for two-sided tests must have tsmethod='central'")
        pval<-switch(alternative,
            less=pbinom(x[1],n,p) - midp.adjustment,
            greater=pbinom(x[1]-1,n,p,lower.tail=FALSE) - midp.adjustment,
            two.sided=exactbinomPvals(x[1],n,p,tsmethod=tsmethod,relErr=relErr,midp=midp)$pvals)
    }
    else {
        m <- r * T
        if (midp){
             if (alternative=="two.sided" & tsmethod!="central") 
                stop("mid-p method for two-sided tests must have tsmethod='central'")
            midp.less<-ppois(x, m)-0.5*dpois(x,m)
            midp.greater<-ppois(x - 1, m, lower.tail = FALSE)-0.5*dpois(x,m)
            pval <- switch(alternative, 
                less = midp.less, 
                greater = midp.greater, 
                two.sided = pmin(rep(1,length(midp.less)),2*midp.less,2*midp.greater))
        } else {
            pval <- switch(alternative, 
                less = ppois(x, m), 
                greater = ppois(x - 1, m, lower.tail = FALSE), 
                two.sided = exactpoissonPvals(x,T,r,relErr, tsmethod=tsmethod)$pvals)
        }
    }
        
    if (newplot){
        if (k==2){
            XLAB<-"Null Hypothesis rate ratio"
        } else {
            XLAB<-"Null Hypothesis rate"
        }
        if (dolog){ LOG<-"x"
        } else LOG<-""
        YLAB<-"p-value"
        if (midp) YLAB<-"mid p-value"
        plot(r,pval,xlab=XLAB,ylab=YLAB,type="n",log=LOG,...)
    } 
    if (dopoints){
        points(r,pval,...)
    } 
    if (dolines){
        lines(r,pval,...)
    }
    if (doci){
        ci<-poisson.exact(x,T,tsmethod=tsmethod,alternative=alternative,
            control=binomControl(relErr=relErr),
            conf.level=conf.level,midp=midp)$conf.int
        alpha<-1-conf.level
        if (alphaline) lines(range(r),c(alpha,alpha),lty=2)
        lines(c(ci[1],ci[1]),c(0,alpha),...)
        lines(c(ci[2],ci[2]),c(0,alpha),...)
    }
    invisible(list(r=r,p.value=pval))
}

#exactpoissonPlot(c(2),c(17877))