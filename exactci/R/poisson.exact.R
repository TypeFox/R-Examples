poisson.exact<-function (x, T = 1, r = 1, alternative = c("two.sided", "less", 
    "greater"),  tsmethod=c("central","minlike","blaker"), conf.level = 0.95,
    control=binomControl(), plot=FALSE, midp=FALSE) 
{
    ## these values should not need to be changed for most uses of the function
    relErr<-control$relErr
    tol<-control$tol
    pRange<-control$pRange

    tsmethod<-match.arg(tsmethod)
    if (tsmethod!="central" & tsmethod!="minlike" & tsmethod!="blaker") stop("tsmethod must be one of 'central', 'minlike', or 'blaker' ")

    ## copy much of this from poisson.test 
    ## for the two sample case primarily just change binom.test to binom.exact


    DNAME <- deparse(substitute(x))
    DNAME <- paste(DNAME, "time base:", deparse(substitute(T)))
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
    if (!missing(r) && (length(r) > 1 || is.na(r) || r < 0)) 
        stop("'r' must be a single positive number")
    alternative <- match.arg(alternative)
    if (k == 2) {
        RVAL <- binom.exact(x, sum(x), r * T[1]/(r * T[1] + T[2]), 
            alternative = alternative, tsmethod=tsmethod, conf.level=conf.level,
            control=control, midp=midp)
        RVAL$data.name <- DNAME
        RVAL$statistic <- x[1]
        RVAL$parameter <- sum(x) * r * T[1]/sum(T * c(1, r))
        names(RVAL$statistic) <- c("count1")
        names(RVAL$parameter) <- c("expected count1")
        RVAL$estimate <- (x[1]/T[1])/(x[2]/T[2])
        names(RVAL$estimate) <- "rate ratio"
        pp <- RVAL$conf.int
        RVAL$conf.int <- pp/(1 - pp) * T[2]/T[1]
        ## Sept 22, 2012, forgot to transform precision also
        if (tsmethod=="blaker" | tsmethod=="minlike"){
            pplowerprec<- attr(RVAL$conf.int,"conf.limit.prec")$lower
            rval.lower.prec<-pplowerprec/(1-pplowerprec)*T[2]/T[1]
            ppupperprec<- attr(RVAL$conf.int,"conf.limit.prec")$upper
            rval.upper.prec<-ppupperprec/(1-ppupperprec)*T[2]/T[1]
            attr(RVAL$conf.int,"conf.limit.prec")<-list(lower=rval.lower.prec,
                upper=rval.upper.prec)
        }
        names(r) <- "rate ratio"
        RVAL$null.value <- r
        methodphrase<-"Exact one-sided Comparison of Poisson rates"
        if (alternative=="two.sided")
            methodphrase<-switch(tsmethod,
                minlike="Exact two-sided Poisson test (sum of minimum likelihood method)",
                central="Exact two-sided Poisson test (central method)",
                blaker="Exact two-sided Poisson test (Blaker's method)")
        if (midp) methodphrase<-paste0(methodphrase,", mid-p version")

        RVAL$method <- methodphrase
    }
    else {
        m <- r * T
        if (midp){
            if (alternative=="two.sided" & tsmethod!="central") stop("midp=TRUE only available for tsmethod='central'")
            midp.less<-ppois(x, m)-0.5*dpois(x,m)
            midp.greater<-ppois(x - 1, m, lower.tail = FALSE)-0.5*dpois(x,m)
            PVAL <- switch(alternative, 
                less = midp.less,
                greater = midp.greater, 
                two.sided =pmin(rep(1,length(midp.less)), 2*midp.less, 2*midp.greater)) 
            p.L<-function(x,alpha){
                if (x==0){
                    out<- 0
                } else  {
                    rootfunc<-function(mu){
                        # check function without midp correction 
                        # against usual pois.exact()$conf.int
                        #1-ppois(x-1,mu) - alpha
                        # with midp correction
                        1-ppois(x,mu)+0.5*dpois(x,mu) - alpha
                    }
                    out<-uniroot(rootfunc,c(0,qgamma(alpha,x)+1))$root
                }
                out
            }
            p.U<-function(x,alpha){
                rootfunc<-function(mu){
                    # check function without midp correction
                    # against usual binom.test()$conf.int
                    #ppois(x,mu) - alpha
                    # with midp correction
                    ppois(x,mu)-0.5*dpois(x,mu) - alpha
                }
                uniroot(rootfunc,
                    c(0,qgamma(1 - alpha, x + 1)+1))$root
            }
        } else {
            PVAL <- switch(alternative, 
                less = ppois(x, m),
                greater = ppois(x - 1, m, lower.tail = FALSE), 
                two.sided = exactpoissonPval(x,T,r,relErr, tsmethod=tsmethod))
            p.L <- function(x, alpha) {
                if (x == 0) 
                    0
                else qgamma(alpha, x)
            }
            p.U <- function(x, alpha) qgamma(1 - alpha, x + 1)
        }

        if (alternative=="less"){
            CINT<-c(0, p.U(x, 1 - conf.level))
        } else if (alternative=="greater"){
            CINT<-c(p.L(x, 1 - conf.level), Inf)
        } else {
            if (tsmethod=="central"){
                alpha <- (1 - conf.level)/2 
                CINT<-c(p.L(x, alpha), p.U(x, alpha))
            } else {
                if (midp) stop("midp=TRUE only available for tsmethod='central'")
                CINT<-exactpoissonCI(x,tsmethod=tsmethod,conf.level=conf.level)
            }
        }
        CINT<-CINT/T
        attr(CINT,"conf.level")<-conf.level

        ESTIMATE <- x/T
        names(x) <- "number of events"
        names(T) <- "time base"
        names(ESTIMATE) <- names(r) <- "event rate"
        methodphrase<-"Exact one-sided Poisson rate test"
        if (alternative=="two.sided")
            methodphrase<-switch(tsmethod,
                minlike="Exact two-sided Poisson test (sum of minimum likelihood method)",
                central="Exact two-sided Poisson test (central method)",
                blaker="Exact two-sided Poisson test (Blaker's method)")

        if (midp) methodphrase<-paste0(methodphrase,", mid-p version")

        RVAL<-structure(list(statistic = x, parameter = T, p.value = PVAL, 
            conf.int = CINT, estimate = ESTIMATE, null.value = r, 
            alternative = alternative, method = methodphrase, 
            data.name = DNAME), class = "htest")
    }

    if (plot){
        ## use all defaults for plot
        lowerRange<-ifelse(RVAL$conf.int[1]==0, .8*min(r,RVAL$conf.int[2]), .8*min(RVAL$conf.int[1],r) )
        upperRange<-ifelse(RVAL$conf.int[2]==Inf, 1.2*max(r,RVAL$conf.int[1]), 1.2*max(RVAL$conf.int[2],r) )
        exactpoissonPlot(x,T,tsmethod=tsmethod,rRange=c(lowerRange,upperRange),alternative=alternative,conf.level=conf.level,col="gray")
        points(r,RVAL$p.value,col="black")
    }
    return(RVAL)
}
#poisson.exact(c(5,4),c(132412,311312))