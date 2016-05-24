###  To fit a smooth curve to histogram-type data using FMKL GLD

## We use percentile method (or quantile matching method to fit FMKL
## GLD to find initial values, then find MLE numerically.

## Created on 2014/06/12 

## a hybrid GLD fitting method:
fit.GLD <- function(x, lbound, ubound, method='chisquare')
{
    res <- NULL; pv0 <- 0;
    out <- fit.GLD.FMKL(x=x,lbound=lbound,ubound=ubound,
                        percentile='appro',mle=FALSE)
    if(any(!is.finite(out$pars)) || out$pars[2] <= 0){
        warning("MOP (approximate) failed")
    }else{
        res <- out
        pv0 <- gof(out, method=method)$p.value
    }
    out <- fit.GLD.FMKL(x=x,lbound=lbound,ubound=ubound,
                        percentile='appro',mle=TRUE)
    if(any(!is.finite(out$pars)) || out$pars[2] <= 0){
        warning("MLE (approximate) failed")
    }else{
        pv <- gof(out, method=method)$p.value
        if(pv > pv0){
            res <- out
            pv0 <- pv
        }
    }
    out <- fit.GLD.FMKL(x=x,lbound=lbound,ubound=ubound,
                        percentile='exact',mle=FALSE)
    if(any(!is.finite(out$pars)) || out$pars[2] <= 0){
        warning("MOP (exact) failed")
    }else{
        pv <- gof(out, method=method)$p.value
        if(pv > pv0){
            res <- out
            pv0 <- pv
        }
    }
    out <- fit.GLD.FMKL(x=x,lbound=lbound,ubound=ubound,
                        percentile='exact',mle=TRUE)
    if(any(!is.finite(out$pars)) || out$pars[2] <= 0){
        warning("MLE (exact) failed")
    }else{
        pv <- gof(out, method=method)$p.value
        if(pv > pv0){
            res <- out
            pv0 <- pv
        }
    }
    if(!missing(lbound) && !missing(ubound)){
        out <- fit.GB(x=x,lbound=lbound,ubound=ubound)
        pv <- gof(out, method=method)$p.value
        if(pv > pv0){
            res <- out
            pv0 <- pv
        }
    }
        
    res
}

fit.GLD.FMKL <- function(x, lbound, ubound, percentile='exact', mle=FALSE)
{
    x.raw <- FALSE
    if(class(x) == 'numeric'){
        y <- x; #store raw data to compute percentiles
        x <- hist(x, plot=FALSE)
        x.raw <- TRUE
    }
    if(class(x) == 'histogram')
        x <- binning(x)
    if(class(x) != 'bdata')
        stop("The current version support histogram-type data only.")

    percentile <- match.arg(tolower(percentile),
                        c("exact","approximate"))
    if(x.raw) percentile <- "raw"

    if(missing(lbound)){
        lbound <- NULL
        a <- x$breaks[1]
    }else if(is.numeric(lbound)){
        if(lbound > x$breaks[1])
            stop("'lbound' larger than the lower limit of the first class")
        if(is.finite(lbound))
            a <- lbound
        else
            a <- x$breaks[1]
    }else stop("'lbound' is not numeric")

    if(missing(ubound)){
        ubound <- NULL
        b <- x$breaks[x$nclass+1]
    }else if(is.numeric(ubound)){
        if(ubound < x$breaks[x$nclass+1])
            stop("'ubound' smaller than the upper limit of the last class")
        if(is.finite(ubound))
            b <- ubound
        else
            b <- x$breaks[x$nclass+1]
    }else stop("'ubound' is not numeric")

    
    if(mle){
        ##  if not initial values can be used from percentile-matching
        ##  fitting, use boundary information to improve the fitting.
        res <- .gld.qtlmatching(x, FALSE,a,b,lbound,ubound)
        if(is.null(res)){
            xsum <- .summary.bdata(x)
            res <- c(xsum$mean,sqrt(xsum$var),1,1)
        }
        res <- .gld.mle(res, x, a, b, lbound, ubound)
    }else{
        if(x.raw){
            qlevels <- c(0.1, 0.25, 0.50, 0.75, 0.90)
            qtls <- as.numeric(quantile(y, qlevels))
            res <- .gld.mop(qlevels,qtls,a,b,lbound,ubound)
        }else{
            x.exact <- switch(percentile, 'exact'=TRUE,'approximate'=FALSE)
            res <- .gld.qtlmatching(x, x.exact,a,b,lbound,ubound)
        }
    }
    
    if(mle)
        method <- "MLE"
    else if(percentile=="exact")
        method <- "Percentile.Matching(exact)"
    else
        method <- "Percentile.Matching(approx)"
    
    structure(
        list(xhist=x, pars=res,method=method,
             call = match.call()),
        class="FMKL")
}

print.FMKL <- function(x,...){
    ## show mean/sd/range/median/skewness/Kurtosis
    cat("FMKL-GLD\n  Call:  ", deparse(x$call), sep = "")
    if(is.null(x$pars)){
        warning("GLD fitting failed!")
    }else{
        cat("\n  Parameters: = (", x$pars[1], ", ",
            x$pars[2], ", ", x$pars[3], ", ",
            x$pars[4],")\n", sep='')
        x.range <- qgld(c(0,0.25,0.5,0.75,1), x$pars) 
        cat("  Range: =[",x.range[1],", ", x.range[5],"]\tMedian: =",
            x.range[3], "\tIQR: =", x.range[4]-x.range[2],
            "\n",sep='')
    }
}

plot.FMKL <- function(x,hist=TRUE,...){
    if(is.null(x$pars))
        warning("GLD fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)
        f0 <- dgld(x0, x$pars)
#        x.sele <- f0 > 0
#        x0 <- x0[x.sele]
#        f0 <- f0[x.sele]
        if(hist){
            plot(res,...)
            lines(f0~x0,...)
        }else{
            plot(f0~x0, ...)
        }
    }
}

lines.FMKL <- function(x,...){
    if(is.null(x$pars))
        warning("GLD fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)
        f0 <- dgld(x0, x$pars)
#        x.sele <- f0 > 0
#        x0 <- x0[x.sele]
#        f0 <- f0[x.sele]
        lines(f0~x0,...)
    }
}

    
.gld.qtlmatching <- function(x, exact,a,b,lbound,ubound){
    a <- x$breaks[1]; b <- x$breaks[x$nclass+1]
    qout <- .summary.bdata(x)
    if(exact){
        qlevels <- qout$exact$levels
        qtls <- qout$exact$qtls
        if(any(is.na(qtls))||length(qtls)<5){
            warning("Exact quantiles not found. GLD fitting failed!")
            lambdas <- NULL
        }else{
            if(!is.null(lbound)){
                if(is.finite(lbound)){
                    qlevels[1] <- 1e-18
                    qtls[1] <- lbound
                }
            }
            if(!is.null(ubound)){
                if(is.finite(ubound)){
                    qlevels[5] <- 1-1e-18
                    qtls[5] <- ubound
                }
            }
            lambdas <- .gld.mop(qlevels,qtls,a,b,lbound,ubound)
        }
    }else{
        qlevels <- qout$approx$levels
        qtls <- qout$approx$qtls
        if(any(is.na(qtls))||length(qtls)<5){
            warning("Approximate quantiles not available. GLD fitting failed!")
            lambdas <- NULL
        }else{
            if(!is.null(lbound)){
                if(is.finite(lbound)){
                    qlevels[1] <- 1e-18
                    qtls[1] <- lbound
                }
            }
            if(!is.null(ubound)){
                if(is.finite(ubound)){
                    qlevels[5] <- 1-1e-18
                    qtls[5] <- ubound
                }
            }
            lambdas <- .gld.mop(qlevels,qtls,a,b,lbound,ubound)
        }
    }
    lambdas
}


.gld.mle <- function(pars, x,a,b,lbound,ubound){
    .g0.mle <- function(lmd){
        breaks <- x$breaks
        counts <- x$counts
        Fx <- diff(pgld(breaks, lmd))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }
    .g1.mle <- function(lmd){
        l1 <- lbound + 1/(lmd[1]*lmd[2])
        lambdas <- c(l1, lmd)
        breaks <- x$breaks
        counts <- x$counts
        Fx <- diff(pgld(breaks, lambdas))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }
    .g2.mle <- function(lmd){
        l1 <- ubound - 1/(lmd[1]*lmd[3])
        lambdas <- c(l1, lmd)
        breaks <- x$breaks
        counts <- x$counts
        Fx <- diff(pgld(breaks, lambdas))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }
    .g3.mle <- function(lmd){
        l3 <- lmd[1]; l4 <- lmd[2]
        l2 <- (1/l3+1/l4)/(ubound-lbound)
        l1 <- ubound - 1/(l2*l4)
        lambdas <- c(l1, l2,l3,l4)
        breaks <- x$breaks
        counts <- x$counts
        Fx <- diff(pgld(breaks, lambdas))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }

    
    ##  constraints:
    if(is.null(lbound)){
        l3l <- -Inf; l3u <- Inf; lfinite <- FALSE
    }else if(is.finite(lbound)){
        l3l <- 1e-8; l3u <- Inf; lfinite <- TRUE
    }else{
        l3l <- -Inf; l3u <- -1e-8; lfinite <- FALSE
    }
    
    if(is.null(ubound)){
        l4l <- -Inf; l4u <- Inf; ufinite <- FALSE
    }else if(is.finite(ubound)){
        l4l <- 1e-8; l4u <- Inf; ufinite <- TRUE
    }else{
        l4l <- -Inf; l4u <- -1e-8; ufinite <- FALSE
    }
    ## find the MLEs
    if(lfinite&&ufinite){
        out <- optim(pars[-c(1,2)], .g3.mle, method="L-BFGS-B",
                     lower=c(l3l, l4l),
                     upper=c(l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
        l3 <- pars[1]; l4 <- pars[2]
        l2 <- (1/l3+1/l4)/(ubound-lbound)
        l1 <- ubound - 1/(l2*l4)
        pars <- c(l1, l2,l3,l4)
    }
    
    if(!lfinite&&ufinite){
        out <- optim(pars[-1], .g2.mle, method="L-BFGS-B",
                     lower=c(1e-10,l3l, l4l),
                     upper=c(Inf,l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
        l1 <- ubound - 1/(pars[1]*pars[3])
        pars <- c(l1, pars)
    }

    if(lfinite&&!ufinite){
        out <- optim(pars[-1], .g1.mle, method="L-BFGS-B",
                     lower=c(1e-10,l3l, l4l),
                     upper=c(Inf,l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
        l1 <- lbound + 1/(pars[1]*pars[2])
        pars <- c(l1, pars)
    }

    if(!lfinite&&!ufinite){
        out <- optim(pars, .g0.mle, method="L-BFGS-B",
                      lower=c(-Inf,1e-10,l3l, l4l),
                      upper=c(Inf,Inf,l3u,l4u),
                     control = list(maxit=10000))
        pars <- as.numeric(out$par)
    }



    if(pars[2] <= 0){
        warning("No valid parameters found")
        pars <- NULL
    }
    pars
}


## l1 can be computed based on the median, which is more stable.  So
## we can just adjust lambda2 (dispersion parameter to cover the whole
## range).

.fexp <- function(p,lmd){
    stopifnot(p>=0&&p<=1)
    if(lmd==0)
        res <- 1e10
    else{
        if(p==0){
            if(lmd>=0)
                res <- 0
            else
                res <- 1e10
        }else
            res <- (p^lmd -1)/lmd
    }
}

#        if(!is.finite(out)){
#            cat("\nlmd:", lmd,
#                "\n a=:",a, "\tb=:",b,
#                "\n p=",p,
#                "\n q=",q,
#                "\nRhos", c(rho3hat, rho4hat, rho3, rho4), "\n")
#            out <- 1e10
#        }

.gld.mop <- function(p,q,a,b,lbound,ubound){
    .g.mop <- function(lmd){
        rho3hat <- (q[3] - q[1])/(q[5] - q[3])
        rho4hat <- (q[4]-q[2])/(q[5]-q[1])
        Q1 <- .fexp(p[1],lmd[1]) - .fexp(1-p[1],lmd[2])
        Q2 <- .fexp(p[2],lmd[1]) - .fexp(1-p[2],lmd[2])
        Q3 <- .fexp(p[3],lmd[1]) - .fexp(1-p[3],lmd[2])
        Q4 <- .fexp(p[4],lmd[1]) - .fexp(1-p[4],lmd[2])
        Q5 <- .fexp(p[5],lmd[1]) - .fexp(1-p[5],lmd[2])

        rho3 <- ifelse(Q5-Q3==0, 1e10, (Q3-Q1)/(Q5-Q3))
        rho4 <- ifelse(Q5-Q1==0, 1e10, (Q4-Q2)/(Q5-Q1))
        out <- (rho3-rho3hat)^2+(rho4-rho4hat)^2
    }
    lambdas <- c(1,1)

    lb <- FALSE
    if(is.null(lbound)){
        l3l <- -Inf; l3u <- Inf; 
    }else if(is.finite(lbound)){
        l3l <- 1e-8; l3u <- Inf; lb <- TRUE
    }else{
        l3l <- -Inf; l3u <- -1e-8;
    }

    ub <- FALSE
    if(is.null(ubound)){
        l4l <- -Inf; l4u <- Inf
    }else if(is.finite(ubound)){
        l4l <- 1e-8; l4u <- Inf; ub <- TRUE
    }else{
        l4l <- -Inf; l4u <- -1e-8
    }

    pars <- optim(lambdas, .g.mop,method="L-BFGS-B",
                  lower=c(l3l,l4l),upper=c(l3u,l4u),
                  control = list(maxit=10000))$par
    
    l3 <- pars[1]; l4 <- pars[2]
    Q5 <- (p[5]^l3-1)/l3 - ((1-p[5])^l4-1)/l4
    Q3 <- (p[3]^l3-1)/l3 - ((1-p[3])^l4-1)/l4
    Q1 <- (p[1]^l3-1)/l3 - ((1-p[1])^l4-1)/l4
    l2 <- (Q5-Q1)/(q[5]-q[1])
    l1 <- q[3] - Q3/l2

    if(lb){
        if(ub){
            l2 <- (1/l3+1/l4)/(ubound-lbound)
            l1 <- ubound - 1/(l2*l4)
        }else{
            l1 <- lbound + 1/(l2*l3)
        }
    }else{
        if(ub)
            l1 <- ubound - 1/(l2*l4)
    }
    res <- c(l1,l2,l3,l4)
#    l20 <- l2
#    ## check whether adjustment is needed
#    xmax <- qgld(1,res)
#    xmin <- qgld(0,res)

#    if(a < xmin){
#        l2a <- 1/((l1-a)*l3)
#        if(l2a > 0) l2 <- l2a
#    }
    
#    if(b > xmax){
#        l2b <- 1/((b-l1)*l4)
#        if(l2b > 0) l2 <- min(l2, l2b)
#    }

#    res <- c(l1,l2,l3,l4)
#    xmax <- qgld(1,res)
#    xmin <- qgld(0,res)
#    if(a < xmin|| b > xmax){
#        l2c <- (1/l4+1/l3)/(b-a)
#        if(l2c > 0) l2 <- min(l2, l2c)
#    }
#    if(l20 != l2)
#        warning("Ranges not match. Estimates adjusted.")

#    res <- c(l1,l2,l3,l4)
#    res
#    res <- c(l1,l20,l3,l4)
    
}

.gld.nls <- function(x){
    if(class(x)=='histogram'){
        x <- binning(x)
    }
    if(class(x)=='bdata'){
        Fn <- cumsum(x$counts)/sum(x$counts)
        Fn <- Fn[-x$nclass]
        xn <- x$breaks[-c(1,x$nclass)]
    }else{#if raw data
        x.out <- .edf(x)
        n <- length(x.out$y)
        Fn <- x.out$y[-n]
        xn <- x.out$x[-n]
    }
    nlmod0 <- nls(xn~l1+((Fn^l3-1)/l3-((1-Fn)^l4-1)/l4)/l2,
                 algorithm='port',
                 lower=c(-Inf, 0.000000001, -Inf,-Inf))

    nlmod <- nls(xn~l1+((Fn^l3-1)/l3-((1-Fn)^l4-1)/l4)/l2,
                 algorithm='port',start=coef(nlmod0),
                 lower=c(-Inf, 0.000000001, -Inf,-Inf))

    pars <- coef(nlmod)
}

rgld <- function(n, lambdas){
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    l1 <- lambdas[1]
    l2 <- lambdas[2]
    l3 <- lambdas[3]
    l4 <- lambdas[4]
    u <- runif(n)
    l1+((u^l3-1)/l3-((1-u)^l4-1)/l4)/l2
}

qgld <- function(p, lambdas){
    if(p < 0 || p > 1)
        stop("Invalid 'p' value")
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    l1 <- lambdas[1]
    l2 <- lambdas[2]
    l3 <- lambdas[3]
    l4 <- lambdas[4]
    u <- p
    l1+((u^l3-1)/l3-((1-u)^l4-1)/l4)/l2
}

.fungld <- function(u,q,lambdas){
    qgld(u,lambdas) - q
}

.gld.proot <- function(q, lambdas){
    if(is.null(lambdas))
        stop("Empty GLD parameters")
    if(is.null(q))
        stop("'q' value missing")
    if(lambdas[2]<=0){
        stop("Invalid GLD parameters")
    }
    xmax <- qgld(1,lambdas)
    xmin <- qgld(0,lambdas)
    if(q == Inf)
        res <- 1
    else if(q==-Inf)
        res <- 0
    else{
        if(q <= xmin)
            res <- 0
        else if(q >= xmax)
            res <- 1
        else{
            res <- .gld.root(q,lambdas=lambdas)
        }
    }
    res
}

pgld <- function(q, lambdas){
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    as.numeric(lapply(q, .gld.proot, lambdas=lambdas))
}


.gld.droot <- function(q, lambdas){
    xmax <- qgld(1,lambdas)
    xmin <- qgld(0,lambdas)
    if(!is.finite(q)){
        res <- 0
    }else if(q < xmin)
        res <- 0
    else if(q > xmax)
        res <- 0
    else{
        u <- .gld.root(q,lambdas=lambdas)
        l1 <- lambdas[1]
        l2 <- lambdas[2]
        l3 <- lambdas[3]
        l4 <- lambdas[4]
        res <- u^(l3-1)+(1-u)^(l4-1)
        res <- l2/res
    }
    res
}

dgld <- function(x, lambdas){
    if(any(!is.finite(lambdas)) || lambdas[2]<=0)
        stop("Invalid FMKL-GLD parameters")
    as.numeric(lapply(x, .gld.droot, lambdas=lambdas))
}

.gld.root <- function(q, lambdas){
    .Fortran(.F_rootGldFmklBisection,
             u=as.double(q), as.double(lambdas))$u
}


