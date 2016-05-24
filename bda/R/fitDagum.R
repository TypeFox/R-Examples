###  Dagum distribution (type I)

## F(x,a,b,p) = (1+(x/b)^(-a))^(-p)

## Created on 2014/06/12 

fit.Dagum <- function(x)
{
    if(class(x)=="numeric")
        x <- hist(x, plot=FALSE)
    if(class(x) == 'histogram')
        x <- binning(x)
    if(class(x) != 'bdata')
        stop("The current version support histogram-type data only.")

    par0 <- .dagum.pm(x)
    res <- .dagumBMLE(x, par0)
    
    structure(
        list(xhist=x, pars=res, 
             call = match.call()),
        class="Dagum")
}

.dagumBMLE <- function(x, pars)
{
    .dagum.bllk <- function(pars){
        res <- -sum(counts * log(diff(pdagum(breaks, pars[1], pars[2], pars[3]))))
        if(!is.finite(res)) res <- 99999999999.
        res
    }

    breaks <- x$breaks
    counts <- x$counts
    
    tol <- 0.0000001
    out <- optim(pars, .dagum.bllk, ## gr=.weibull.grbllk,
                 method="L-BFGS-B", lower=c(tol,tol, tol))
    as.numeric(out$par)
}

.dagum.pm <- function(x){
    ## to estimate the parameters of Dagum using percentile matching
    ## method.

    .fdagum <- function(p){
        l3 <- log(exp(-log(p3)/p)-1)
        l2 <- log(exp(-log(p2)/p)-1)
        l1 <- log(exp(-log(p1)/p)-1)
        d32 <- log(q3) - log(q2)
        d21 <- log(q2) - log(q1)
        abs((l3-l2)/(l2-l1) - d32/d21)
    }

    y <- .summary.bdata(x)
    if(y$exact){
        k <- 5
        q1 <- y$qtls2[2];
        q2 <- y$qtls2[3]
        q3 <- y$qtls2[4]
        p1 <- y$qlevels2[2];
        p2 <- y$qlevels2[3];
        p3 <- y$qlevels2[4];
    }else{
        k <- length(y$qtls1)
        if(k==5){
            q1 <- y$qtls1[2];
            q2 <- y$qtls1[3]
            q3 <- y$qtls1[4]
            p1 <- y$qlevels1[2];
            p2 <- y$qlevels1[3];
            p3 <- y$qlevels1[4];
        }else if(k==4||k==3){
            q1 <- y$qtls1[1];
            q2 <- y$qtls1[2]
            q3 <- y$qtls1[3]
            p1 <- y$qlevels1[1];
            p2 <- y$qlevels1[2];
            p3 <- y$qlevels1[3];
        }
    }
    if(k < 3){
        res <- c(1, x$breaks[x$nclass]-x$breaks[2], 1)
    }else{
        pars <- 1
        out <- optim(pars, .fdagum, method="L-BFGS-B",lower=1.e-8)
        p <- as.numeric(out$par)

        l3 <- log(exp(-log(p3)/p)-1)
        l2 <- log(exp(-log(p2)/p)-1)
        l1 <- log(exp(-log(p1)/p)-1)
        d32 <- log(q3) - log(q2)
        a <- -(l3-l2)/d32
        b <- exp(l3/a + log(q3))
        res <- c(a,b,p)
    }
    res
}

print.Dagum <- function(x,...){
    cat("The Dagum Distribution (type I)\n  Call:  ", deparse(x$call), sep = "")
    if(is.null(x$pars)){
        warning("Dagum fitting failed!")
    }else{
        a <- x$pars[1]
        b <- x$pars[2]
        p <- x$pars[3]
        cat("\n  Parameters: = (", a, ", ",
            b, ", ", p,")\n", sep='')
        cat("  Range: =[ 0, Inf)\tMedian: =",
            qdagum(0.5, a, b, p), "\n",sep='')
    }
}

plot.Dagum <- function(x,hist=TRUE,...){
    if(is.null(x$pars))
        warning("GLD fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)
        a <- x$pars[1]
        b <- x$pars[2]
        p <- x$pars[3]
        f0 <- ddagum(x0, a,b,p)
        if(hist){
            plot(res,...)
            lines(f0~x0,...)
        }else{
            plot(f0~x0, ...)
        }
    }
}

lines.Dagum <- function(x,...){
    if(is.null(x$pars))
        warning("GLD fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)
        a <- x$pars[1]
        b <- x$pars[2]
        p <- x$pars[3]
        f0 <- ddagum(x0, a,b,p)
        lines(f0~x0,...)
    }
}

pdagum <- function(q, shape1,scale,shape2){
    stopifnot(shape1>0&&scale>0&shape2>0)
    res <- rep(0, length(q))
    x.sele1 <-  q <= 0
    if(any(x.sele1))
        res[x.sele1] <- 0
    x.sele2 <- !is.finite(q) & q > 0  
    if(any(x.sele2))
        res[x.sele2] <- 1
    sele <- !x.sele1 & !x.sele2
    if(any(sele)){
        q <- q[sele]
        a <- shape1
        b <- scale
        p <- shape2
        out <- (1+(q/b)^(-a))^(-p)
        res[sele] <- out 
    }
    res
}

ddagum <- function(x, shape1,scale,shape2){
    stopifnot(shape1>0&&scale>0&shape2>0)
    res <- rep(0, length(x))
    x.sele <-  x <= 0 | !is.finite(x)
    if(any(x.sele))
        res[x.sele] <- 0
    if(any(!x.sele)){
        x <- x[!x.sele]
        a <- shape1
        b <- scale
        p <- shape2
        xb <- x/b
        out <- a*p/x*(xb)^(a*p)*(1+(xb)^a)^(-1-p)
        res[!x.sele] <- out
    }
    res
}

qdagum <- function(p, shape1,scale,shape2){
    stopifnot(shape1>0&&scale>0&shape2>0)
    res <- rep(0, length(p))
    sele1 <- p < 0
    if(any(sele1))
        res[sele1] <- NA
    sele2 <- p == 0
    if(any(sele2))
        res[sele2] <- 0
    sele3 <- p == 1
    if(any(sele3))
        res[sele3] <- Inf
    sele4 <- p > 1
    if(any(sele4))
        res[sele4] <- NA
    sele <- !sele1 & !sele2 & !sele3 & !sele4
    
    if(any(sele)){
        u <- p[sele]
        a <- shape1
        b <- scale
        p <- shape2
        out <- b*(u^(-1/p)-1)^(-1/a)
        res[sele] <- out
    }
    res
}
